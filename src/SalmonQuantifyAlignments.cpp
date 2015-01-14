
extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
}

// for cpp-format
#include "format.h"

// are these used?
#include <boost/dynamic_bitset.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/lockfree/queue.hpp>

#include <tbb/atomic.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <random>
#include <memory>
#include <exception>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <memory>
#include <condition_variable>

#include <tbb/concurrent_queue.h>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/program_options.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/normal.hpp>

#include "ClusterForest.hpp"
#include "AlignmentLibrary.hpp"
#include "MiniBatchInfo.hpp"
#include "BAMQueue.hpp"
#include "SailfishMath.hpp"
#include "FASTAParser.hpp"
#include "LibraryFormat.hpp"
#include "Transcript.hpp"
#include "ReadPair.hpp"
#include "ErrorModel.hpp"
#include "FragmentLengthDistribution.hpp"
#include "TranscriptCluster.hpp"
#include "SailfishUtils.hpp"
#include "SalmonUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonOpts.hpp"
#include "NullFragmentFilter.hpp"
#include "Sampler.hpp"
#include "spdlog/spdlog.h"

namespace bfs = boost::filesystem;
using sailfish::math::LOG_0;
using sailfish::math::LOG_1;
using sailfish::math::logAdd;
using sailfish::math::logSub;

template <typename FragT>
using AlignmentBatch = std::vector<FragT>;

template <typename FragT>
using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

using PriorAbundanceVector = std::vector<double>;
using PosteriorAbundanceVector = std::vector<double>;

struct RefSeq {
    RefSeq(char* name, uint32_t len) : RefName(name), RefLength(len) {}
    std::string RefName;
    uint32_t RefLength;
};

/**
 * Multiple each element in the vector `vec` by the factor `scale`.
 */
template <typename T>
void scaleBy(std::vector<T>& vec, T scale) {
    std::for_each(vec.begin(), vec.end(), [scale](T& ele)->void { ele *= scale; });
}

template <typename FragT>
void processMiniBatch(AlignmentLibrary<FragT>& alnLib,
                      MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                      std::condition_variable& workAvailable,
                      std::mutex& cvmutex,
                      std::atomic<bool>& doneParsing,
                      std::atomic<size_t>& activeBatches,
                      const SalmonOpts& salmonOpts,
                      bool& burnedIn,
                      bool initialRound,
                      std::atomic<size_t>& processedReads) {

    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

    using sailfish::math::LOG_0;
    using sailfish::math::logAdd;
    using sailfish::math::logSub;

    auto& refs = alnLib.transcripts();
    auto& clusterForest = alnLib.clusterForest();
    auto& fragmentQueue = alnLib.fragmentQueue();
    auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();
    auto& fragLengthDist = alnLib.fragmentLengthDistribution();
    auto& errMod = alnLib.errorModel();

    std::chrono::microseconds sleepTime(1);
    MiniBatchInfo<AlignmentGroup<FragT*>>* miniBatch;
    bool updateCounts = initialRound;
    size_t numTranscripts = refs.size();

    while (!doneParsing) {
        miniBatch = nullptr;
        {
            std::unique_lock<std::mutex> l(cvmutex);
            workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() { return workQueue.try_pop(miniBatch) or doneParsing; });
        }
        if (miniBatch != nullptr) {
            ++activeBatches;
            size_t batchReads{0};
            double logForgettingMass = miniBatch->logForgettingMass;
            std::vector<AlignmentGroup<FragT*>*>& alignmentGroups = *(miniBatch->alignments);

            using TranscriptID = size_t;
            using HitIDVector = std::vector<size_t>;
            using HitProbVector = std::vector<double>;

            std::unordered_map<TranscriptID, std::vector<FragT*>> hitList;
            for (auto alnGroup : alignmentGroups) {
                for (auto a : alnGroup->alignments()) {
                    auto transcriptID = a->transcriptID();
                    if (transcriptID < 0 or transcriptID >= refs.size()) {
                        std::cerr << "Invalid Transcript ID: " << transcriptID << "\n";
                    }
                    hitList[transcriptID].emplace_back(a);
                }
            }

            {
                // E-step
                //boost::timer::auto_cpu_timer t(3, "E-step part 2 took %w sec.\n");

                // Iterate over each group of alignments (a group consists of all alignments reported
                // for a single read).  Distribute the read's mass proportionally dependent on the
                // current
                for (auto alnGroup : alignmentGroups) {
                    double sumOfAlignProbs{LOG_0};
                    // update the cluster-level properties
                    bool transcriptUnique{true};
                    auto firstTranscriptID = alnGroup->alignments().front()->transcriptID();
                    std::unordered_set<size_t> observedTranscripts;
                    for (auto& aln : alnGroup->alignments()) {
                        auto transcriptID = aln->transcriptID();
                        auto& transcript = refs[transcriptID];
                        transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                        double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;

                        double logFragProb = sailfish::math::LOG_1;
                        //double fragLength = aln.fragLen();

                        if (!salmonOpts.noFragLengthDist) {
                            if(aln->fragLen() == 0) {
                                if (aln->isLeft() and transcript.RefLength - aln->left() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(transcript.RefLength - aln->left());
                                } else if (aln->isRight() and aln->right() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(aln->right());
                                }
                            } else {
                                logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln->fragLen()));
                            }

                        }

                        // @TODO: handle this case better
                        //double fragProb = cdf(fragLengthDist, fragLength + 0.5) - cdf(fragLengthDist, fragLength - 0.5);
                        //fragProb = std::max(fragProb, 1e-3);
                        //fragProb /= cdf(fragLengthDist, refLength);

                        // Some aligners (e.g. Bowtie) don't output alignment qualities, so prepare for this.
                        /*
                        auto q1 = aln.read1->core.qual;
                        auto q2 = aln.read2->core.qual;
                        double logP1 = (q1 == 255) ? calcQuality(aln.read1) : std::log(std::pow(10.0, -q1 * 0.1));
                        double logP2 = (q2 == 255) ? calcQuality(aln.read2) : std::log(std::pow(10.0, -q2 * 0.1));
                        */

                        // The alignment probability is the product of a transcript-level term (based on abundance and) an alignment-level
                        // term below which is P(Q_1) * P(Q_2) * P(F | T)
                        double logRefLength = std::log(refLength);
                        // double adjustedTranscriptLength = std::max(refLength - aln.fragLen() + 1, 1.0);
                        // double logStartPosProb = std::log(1.0 / adjustedTranscriptLength);

                        // P(Fn | Tn) = Probability of selecting a fragment of this length, given the transcript is t
                        // d(Fn) / sum_x = 1^{lt} d(x)
                        //double logConditionalFragLengthProb = logFragProb  - fragLengthDist.cmf(refLength);
                        //double logProbStartPos = logStartPosProb + logConditionalFragLengthProb;
                        //double qualProb = logProbStartPos + aln.logQualProb();
                        double qualProb = -logRefLength + logFragProb + aln->logQualProb();
                        double transcriptLogCount = transcript.mass(initialRound);

                        if ( transcriptLogCount != LOG_0 ) {

                            double errLike = sailfish::math::LOG_1;

                            if (burnedIn and salmonOpts.useErrorModel) {
                                errLike = errMod.logLikelihood(*aln, transcript);
                            }
                            //aln.logProb = (transcriptLogCount - logRefLength) + qualProb + errLike;
                            aln->logProb = transcriptLogCount + qualProb + errLike;

                            sumOfAlignProbs = logAdd(sumOfAlignProbs, aln->logProb);
                            if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                                refs[transcriptID].addTotalCount(1);
                                observedTranscripts.insert(transcriptID);
                            }
                        } else {
                            aln->logProb = LOG_0;
                        }
                    }
                    // normalize the hits
                    if (sumOfAlignProbs == LOG_0) { std::cerr << "0 probability fragment; skipping\n"; continue; }
                    for (auto& aln : alnGroup->alignments()) {
                        aln->logProb -= sumOfAlignProbs;

                        double r = uni(eng);
                        if (!burnedIn and r < std::exp(aln->logProb)) {
                            auto transcriptID = aln->transcriptID();
                            auto& transcript = refs[transcriptID];
                            if (salmonOpts.useErrorModel) {
                                errMod.update(*aln, transcript, aln->logProb, logForgettingMass);
                            }
                            if (aln->isPaired()) {
                                double fragLength = aln->fragLen();//std::abs(aln.read1->core.pos - aln.read2->core.pos) + aln.read2->core.l_qseq;
                                fragLengthDist.addVal(fragLength, logForgettingMass);
                            }
                        }
                    }
                    // update the single target transcript
                    if (transcriptUnique) {
                        if (updateCounts) {
                            refs[firstTranscriptID].addUniqueCount(1);
                        }
                        clusterForest.updateCluster(firstTranscriptID, 1,
                                                    logForgettingMass, updateCounts);
                    } else { // or the appropriate clusters
                        // ughh . . . C++ still has some very rough edges
                        clusterForest.template mergeClusters<FragT>(alnGroup->alignments().begin(),
                                                           alnGroup->alignments().end());
                        clusterForest.updateCluster(alnGroup->alignments().front()->transcriptID(),
                                                    1, logForgettingMass, updateCounts);
                    }

                    ++batchReads;
                } // end read group
            }// end timer

            double individualTotal = LOG_0;
            {
                // M-step
                for (auto kv = hitList.begin(); kv != hitList.end(); ++kv) {
                    auto transcriptID = kv->first;
                    // The target must be a valid transcript
                    if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr << "index " << transcriptID << " out of bounds\n"; }

                    auto& transcript = refs[transcriptID];

                    // The prior probability
                    double hitMass{LOG_0};

                    // The set of alignments that match transcriptID
                    auto& hits = kv->second;
                    std::for_each(hits.begin(), hits.end(), [&](FragT* aln) -> void {
                            if (!std::isfinite(aln->logProb)) { std::cerr << "hitMass = " << aln->logProb << "\n"; }
                            hitMass = logAdd(hitMass, aln->logProb);
                    });

                    double updateMass = logForgettingMass + hitMass;
                    individualTotal = logAdd(individualTotal, updateMass);

                    // Lock the target
                    transcript.addMass(updateMass);

                   // unlock the target
                } // end for
            } // end timer

            miniBatch->release(fragmentQueue, alignmentGroupQueue);
            delete miniBatch;
            --activeBatches;
            processedReads += batchReads;
            if (processedReads >= 5000000 and !burnedIn) { burnedIn = true; }
        }
    } // nothing left to process
}

/**
  *  Quantify the targets given in the file `transcriptFile` using the
  *  alignments given in the file `alignmentFile`, and write the results
  *  to the file `outputFile`.  The reads are assumed to be in the format
  *  specified by `libFmt`.
  *
  */
template <typename FragT>
bool quantifyLibrary(
        AlignmentLibrary<FragT>& alnLib,
        size_t numRequiredFragments,
        uint32_t numQuantThreads,
        const SalmonOpts& salmonOpts) {

    bool burnedIn{false};

    auto& refs = alnLib.transcripts();
    size_t numTranscripts = refs.size();
    size_t miniBatchSize{1000};
    size_t numObservedFragments{0};

    MiniBatchQueue<AlignmentGroup<FragT*>> workQueue;
    double logForgettingMass{std::log(1.0)};
    double forgettingFactor{0.60};
    size_t batchNum{0};
    bool initialRound{true};

    NullFragmentFilter<FragT>* nff = nullptr;

    // Give ourselves some space
    fmt::print(stderr, "\n\n\n\n");

    while (numObservedFragments < numRequiredFragments) {
        if (!initialRound) {
            if (!alnLib.reset(true, nff)) {
                fmt::print(stderr,
                  "\n\n======== WARNING ========\n"
                  "A provided alignment file "
                  "is not a regular file and therefore can't be read from "
                  "more than once.\n\n"
                  "We observed only {} fragments when we wanted at least {}.\n\n"
                  "Please consider re-running Salmon with these alignments "
                  "as a regular file!\n"
                  "==========================\n\n",
                  numObservedFragments, numRequiredFragments);
                break;
            }
        }
        std::atomic<bool> doneParsing{false};
        std::condition_variable workAvailable;
        std::mutex cvmutex;
        std::vector<std::thread> workers;
        std::atomic<size_t> activeBatches{0};
        std::atomic<size_t> processedReads{0};
        for (uint32_t i = 0; i < numQuantThreads; ++i) {
            workers.emplace_back(processMiniBatch<FragT>,
                    std::ref(alnLib),
                    std::ref(workQueue),
                    std::ref(workAvailable), std::ref(cvmutex),
                    std::ref(doneParsing), std::ref(activeBatches),
                    std::ref(salmonOpts),
                    std::ref(burnedIn),
                    initialRound,
                    std::ref(processedReads));
        }

        size_t numProc{0};

        BAMQueue<FragT>& bq = alnLib.getAlignmentGroupQueue();
        std::vector<AlignmentGroup<FragT*>*>* alignments = new std::vector<AlignmentGroup<FragT*>*>;
        alignments->reserve(miniBatchSize);
        AlignmentGroup<FragT*>* ag;

        bool alignmentGroupsRemain = bq.getAlignmentGroup(ag);
        while (alignmentGroupsRemain or alignments->size() > 0) {
            if (alignmentGroupsRemain) { alignments->push_back(ag); }
            // If this minibatch has reached the size limit, or we have nothing
            // left to fill it up with
            if (alignments->size() >= miniBatchSize or !alignmentGroupsRemain) {
                ++batchNum;
                if (batchNum > 1) {
                    logForgettingMass += forgettingFactor * std::log(static_cast<double>(batchNum-1)) -
                        std::log(std::pow(static_cast<double>(batchNum), forgettingFactor) - 1);
                }
                MiniBatchInfo<AlignmentGroup<FragT*>>* mbi =
                    new MiniBatchInfo<AlignmentGroup<FragT*>>(batchNum, alignments, logForgettingMass);
                workQueue.push(mbi);
                {
                    std::unique_lock<std::mutex> l(cvmutex);
                    workAvailable.notify_one();
                }
                alignments = new std::vector<AlignmentGroup<FragT*>*>;
                alignments->reserve(miniBatchSize);
            }

            if ((numProc % 1000000 == 0) or !alignmentGroupsRemain) {
                fmt::print(stderr, "\r\r{}processed{} {} {}reads in current round{}",
                           ioutils::SET_GREEN, ioutils::SET_RED, numProc,
                           ioutils::SET_GREEN, ioutils::RESET_COLOR);
            }

            ++numProc;
            alignmentGroupsRemain = bq.getAlignmentGroup(ag);
        }
        fmt::print(stderr, "\n");

        // Free the alignments and the vector holding them
        for (auto& aln : *alignments) {
            aln->alignments().clear();
            delete aln; aln = nullptr;
        }
        delete alignments;

        doneParsing = true;
        size_t tnum{0};
        for (auto& t : workers) {
            fmt::print(stderr, "\r\rkilling thread {} . . . ", tnum++);
            {
                std::unique_lock<std::mutex> l(cvmutex);
                workAvailable.notify_all();
            }
            t.join();
            fmt::print(stderr, "done");
        }
        fmt::print(stderr, "\n\n");

        initialRound = false;
        numObservedFragments += alnLib.numMappedReads();

        fmt::print(stderr, "# observed = {} / # required = {}\033[A\033[A\033[A\033[A\033[A",
                   numObservedFragments, numRequiredFragments);
    }

    fmt::print(stderr, "\n\n\n\n");
    return burnedIn;
    // Write the inferred fragment length distribution
    /*
    bfs::path distFileName = outputFile.parent_path() / "flenDist.txt";
    {
        std::unique_ptr<std::FILE, int (*)(std::FILE *)> distOut(std::fopen(distFileName.c_str(), "w"), std::fclose);
        fmt::print(distOut.get(), "{}\n", fragLengthDist.toString());
    }
    */
}

int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    boost::filesystem::path outFilePath,
    bool useStreamingParser,
    size_t numThreads);

int performBiasCorrectionSalmon(
        boost::filesystem::path featureFile,
        boost::filesystem::path expressionFile,
        boost::filesystem::path outputFile,
        size_t numThreads);

int salmonAlignmentQuantify(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    SalmonOpts sopt;

    bool sampleOutput{false};
    bool sampleUnaligned{false};
    bool biasCorrect{false};
    uint32_t numThreads{4};
    size_t requiredObservations{50000000};

    po::options_description generic("salmon quant options");
    generic.add_options()
    ("version,v", "print version string.")
    ("help,h", "produce help message.")
    ("libtype,l", po::value<std::string>()->required(), "Format string describing the library type.")
    ("alignments,a", po::value<vector<string>>()->multitoken()->required(), "input alignment (BAM) file(s).")
    //("alignments,a", po::value<string>()->required(), "input alignment (BAM) file(s)")
    ("maxReadOcc,w", po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(200), "Reads \"mapping\" to more than this many places won't be considered.")
    ("targets,t", po::value<std::string>()->required(), "FASTA format file containing target transcripts.")
    ("threads,p", po::value<uint32_t>(&numThreads)->default_value(6), "The number of threads to use concurrently. "
                                            "The alignment-based quantification mode of salmon is usually I/O bound "
                                            "so until there is a faster multi-threaded SAM/BAM parser to feed the "
                                            "quantification threads, one should not expect much of a speed-up beyond "
                                            "~6 threads.")
    ("useReadCompat,e", po::bool_switch(&(sopt.useReadCompat))->default_value(false), "[Currently Experimental] : "
                        "Use the orientation in which fragments were \"mapped\"  to assign them a probability.  For "
                        "example, fragments with an incorrect relative oritenation with respect  to the provided library "
                        "format string will be assigned a 0 probability.")
    ("noEffectiveLengthCorrection", po::bool_switch(&(sopt.noEffectiveLengthCorrection))->default_value(false), "Disables "
                        "effective length correction when computing the probability that a fragment was generated "
                        "from a transcript.  If this flag is passed in, the fragment length distribution is not taken "
                        "into account when computing this probability.")
    ("noFragLengthDist", po::bool_switch(&(sopt.noFragLengthDist))->default_value(false), "[Currently Experimental] : "
                        "Don't consider concordance with the learned fragment length distribution when trying to determine "
                        "the probability that a fragment has originated from a specified location.  Normally, Fragments with "
                         "unlikely lengths will be assigned a smaller relative probability than those with more likely "
                        "lengths.  When this flag is passed in, the observed fragment length has no effect on that fragment's "
                        "a priori probability.")
    ("useErrorModel", po::bool_switch(&(sopt.useErrorModel))->default_value(false), "[Currently Experimental] : "
                        "Learn and apply an error model for the aligned reads.  This takes into account the "
                        "the observed frequency of different types of mismatches when computing the likelihood of "
                        "a given alignment.")
    ("maxReadLen,r", po::value<uint32_t>(&(sopt.maxExpectedReadLen))->default_value(250), "The maximum expected length of an observed read.  "
                        "This is used to determine the size of the error model.  If a read longer than this is "
                        "encountered, then the error model will be disabled")
    ("output,o", po::value<std::string>()->required(), "Output quantification directory.")
    ("sampleOut,s", po::bool_switch(&sampleOutput)->default_value(false), "Write a \"postSample.bam\" file in the output directory "
                        "that will sample the input alignments according to the estimated transcript abundances. If you're "
                        "going to perform downstream analysis of the alignments with tools which don't, themselves, take "
                        "fragment assignment ambiguity into account, you should use this output.")
    ("sampleUnaligned,u", po::bool_switch(&sampleUnaligned)->default_value(false), "In addition to sampling the aligned reads, also write "
                        "the un-aligned reads to \"posSample.bam\".")
    ("biasCorrect", po::value(&biasCorrect)->zero_tokens(), "[Experimental]: Output both bias-corrected and non-bias-corrected "
                                                             "qunatification estimates.")
    ("numRequiredObs,n", po::value(&requiredObservations)->default_value(50000000),
                                        "The minimum number of observations (mapped reads) that must be observed before "
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the "
                                        "input file, then it will be read through multiple times.")
    ("geneMap,g", po::value<std::string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
                                        "Sailfish will output both quant.sf and quant.genes.sf files, where the latter "
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format "
                                        "where each line contains the name of a transcript and the gene to which it belongs "
                                        "separated by a tab.  The extension of the file is used to determine how the file "
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
                                        "format; files with any other extension are assumed to be in the simple format.");

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(generic).run();

        po::store(orderedOptions, vm);

        if (vm.count("help")) {
            std::cout << "Salmon quant (alignment-based)\n";
            std::cout << generic << std::endl;
            std::exit(0);
        }
        po::notify(vm);

        if (numThreads < 2) {
            fmt::print(stderr, "salmon requires at least 2 threads --- "
                               "setting # of threads = 2");
            numThreads = 2;
        }
        sopt.numThreads = numThreads;

        std::stringstream commentStream;
        commentStream << "# salmon (alignment-based) v" << salmon::version << "\n";
        commentStream << "# [ program ] => salmon \n";
        commentStream << "# [ command ] => quant \n";
        for (auto& opt : orderedOptions.options) {
            commentStream << "# [ " << opt.string_key << " ] => {";
            for (auto& val : opt.value) {
                commentStream << " " << val;
            }
            commentStream << " }\n";
        }
        std::string commentString = commentStream.str();
        fmt::print(stderr, "{}", commentString);

        // Verify the geneMap before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("geneMap")) {
            // Make sure the provided file exists
            geneMapPath = vm["geneMap"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                fmt::print(stderr, "Could not find transcript <=> gene map file {} \n"
                           "Exiting now; please either omit the \'geneMap\' option or "
                           "provide a valid file\n", geneMapPath);
                std::exit(1);
            }
        }

        vector<string> alignmentFileNames = vm["alignments"].as<vector<string>>();
        vector<bfs::path> alignmentFiles;
        for (auto& alignmentFileName : alignmentFileNames) {
            bfs::path alignmentFile(alignmentFileName);//vm["alignments"].as<std::string>());
            if (!bfs::exists(alignmentFile)) {
                std::stringstream ss;
                ss << "The provided alignment file: " << alignmentFile <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            } else {
                alignmentFiles.push_back(alignmentFile);
            }
        }

        std::string libFmtStr = vm["libtype"].as<std::string>();
        LibraryFormat libFmt = sailfish::utils::parseLibraryFormatStringNew(libFmtStr);
        if (libFmt.check()) {
            std::cerr << libFmt << "\n";
        } else {
            std::stringstream ss;
            ss << libFmt << " is invalid!";
            throw std::invalid_argument(ss.str());
        }

        bfs::path outputDirectory(vm["output"].as<std::string>());
        // If the path exists
        if (bfs::exists(outputDirectory)) {
            // If it is not a directory, then complain
            if (!bfs::is_directory(outputDirectory)) {
                std::stringstream errstr;
                errstr << "Path [" << outputDirectory << "] already exists "
                       << "and is not a directory.\n"
                       << "Please either remove this file or choose another "
                       << "output path.\n";
                throw std::invalid_argument(errstr.str());
            }
        } else { // If the path doesn't exist, then create it
            if (!bfs::create_directories(outputDirectory)) { // creation failed for some reason
                std::stringstream errstr;
                errstr << "Could not create output directory ["
                       << outputDirectory << "], please check that it is valid.";
                throw std::invalid_argument(errstr.str());
            }
        }

        bfs::path logDirectory = outputDirectory / "logs";

        // Create the logger and the logging directory
        bfs::create_directory(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "Logs will be written to " << logDirectory.string() << "\n";

        bfs::path logPath = logDirectory / "salmon.log";
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("consoleLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        //g2LogWorker logger(argv[0], logDirectory.string());
        //g2::initializeLogging(&logger);


        if (!sampleOutput and sampleUnaligned) {
            fmt::MemoryWriter wstr;
            wstr << "WARNING: you passed in the (-u/--sampleUnaligned) flag, but did not request a sampled "
                 << "output file (-s/--sampleOut).  This flag will be ignored!\n";
            jointLog->warn() << wstr.str();
        }

        // If we made it this far, the output directory exists
        bfs::path outputFile = outputDirectory / "quant.sf";
        // Now create a subdirectory for any parameters of interest
        bfs::path paramsDir = outputDirectory / "libParams";
        if (!boost::filesystem::exists(paramsDir)) {
            if (!boost::filesystem::create_directory(paramsDir)) {
                fmt::print(stderr, "{}ERROR{}: Could not create "
                           "output directory for experimental parameter "
                           "estimates [{}]. exiting.", ioutils::SET_RED,
                           ioutils::RESET_COLOR, paramsDir);
                std::exit(-1);
            }
        }

        // The transcript file contains the target sequences
        bfs::path transcriptFile(vm["targets"].as<std::string>());

        // Currently, one thread is used for parsing the alignment file.
        // Hopefully, in the future, samtools will implemented multi-threaded
        // BAM/SAM parsing, as this is the current bottleneck.  For the time
        // being, however, the number of quantification threads is the
        // total number of threads - 1.
        uint32_t numParseThreads = std::min(uint32_t(4),
                                            std::max(uint32_t(2), uint32_t(std::ceil(numThreads/2.0))));
        numThreads = std::max(numThreads, numParseThreads);
        uint32_t numQuantThreads = std::max(uint32_t(2), uint32_t(numThreads - numParseThreads));
        sopt.numQuantThreads = numQuantThreads;
        sopt.numParseThreads = numParseThreads;
        std::cerr << "numQuantThreads = " << numQuantThreads << "\n";

        switch (libFmt.type) {
            case ReadType::SINGLE_END:
                {
                    AlignmentLibrary<UnpairedRead> alnLib(alignmentFiles,
                                                          transcriptFile,
                                                          libFmt,
                                                          sopt);
                    bool burnedIn = quantifyLibrary<UnpairedRead>(alnLib, requiredObservations, numQuantThreads, sopt);
                    fmt::print(stderr, "\n\nwriting output \n");
                    salmon::utils::writeAbundances(sopt, alnLib, outputFile, commentString);
                    if (sampleOutput) {
                        bfs::path sampleFilePath = outputDirectory / "postSample.bam";
                        salmon::sampler::sampleLibrary<UnpairedRead>(alnLib, numQuantThreads, sopt, burnedIn, sampleFilePath, sampleUnaligned);
                    }
                }
                break;
            case ReadType::PAIRED_END:
                {
                    AlignmentLibrary<ReadPair> alnLib(alignmentFiles,
                                                      transcriptFile,
                                                      libFmt,
                                                      sopt);
                    bool burnedIn = quantifyLibrary<ReadPair>(alnLib, requiredObservations, numQuantThreads, sopt);
                    fmt::print(stderr, "\n\nwriting output \n");
                    salmon::utils::writeAbundances(sopt, alnLib, outputFile, commentString);

                    // Test writing out the fragment length distribution
                    if (!sopt.noFragLengthDist) {
                        bfs::path distFileName = paramsDir / "flenDist.txt";
                        {
                            std::unique_ptr<std::FILE, int (*)(std::FILE *)> distOut(std::fopen(distFileName.c_str(), "w"), std::fclose);
                            fmt::print(distOut.get(), "{}\n", alnLib.fragmentLengthDistribution().toString());
                        }
                    }

                    if (sampleOutput) {
                        bfs::path sampleFilePath = outputDirectory / "postSample.bam";
                        salmon::sampler::sampleLibrary<ReadPair>(alnLib, numQuantThreads, sopt, burnedIn, sampleFilePath, sampleUnaligned);
                    }
                }
                break;
            default:
                std::cerr << "Cannot quantify library of unknown format "
                          << libFmt << "\n";
                std::exit(1);
        }

        bfs::path estFilePath = outputDirectory / "quant.sf";

        if (biasCorrect) {
            // First, compute the transcript features in case the user
            // ever wants to bias-correct his / her results
            bfs::path transcriptBiasFile(outputDirectory); transcriptBiasFile /= "bias_feats.txt";

            bool useStreamingParser{true};
            std::vector<std::string> transcriptFiles{transcriptFile.string()};
            std::cerr << "computeBiasFeatures( {";
            for (auto& tf : transcriptFiles) {
                std::cerr << "[" << tf << "] ";
            }
            std::cerr << ", " << transcriptBiasFile << ", " << useStreamingParser << ", " << numThreads << ")\n";
            computeBiasFeatures(transcriptFiles, transcriptBiasFile, useStreamingParser, numThreads);

            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrectionSalmon(transcriptBiasFile, estFilePath, biasCorrectedFile, numThreads);
        }

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("geneMap")) {
            try {
                sailfish::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory,
                                                            biasCorrect);
            } catch (std::exception& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                                   "estimates. The gene-level file(s) may not exist",
                                   e.what());
            }
        }

    } catch (po::error& e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "============\n";
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << "============\n";
        std::cerr << argv[0] << " alignment-quant was invoked improperly.\n";
        std::cerr << "For usage information, " <<
            "try " << argv[0] << " quant --help-alignments\nExiting.\n";
        std::exit(1);
    }
    return 0;
}
