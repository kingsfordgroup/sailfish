extern "C" {
#include "htslib/sam.h"
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
void processMiniBatch(MiniBatchQueue<AlignmentGroup<FragT>>& workQueue,
                      std::condition_variable& workAvailable,
                      std::mutex& cvmutex,
                      std::vector<Transcript>& refs,
                      ClusterForest& clusterForest, std::atomic<bool>& doneParsing,
                      std::atomic<size_t>& activeBatches,
                      tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue,
                      tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*>& alignmentGroupQueue,
                      ErrorModel& errMod,
                      FragmentLengthDistribution& fragLengthDist,
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

    std::chrono::microseconds sleepTime(1);
    MiniBatchInfo<AlignmentGroup<FragT>>* miniBatch;
    bool updateCounts = initialRound;
    size_t numTranscripts = refs.size();

    while (!doneParsing) {
        miniBatch = nullptr;
        {
            std::unique_lock<std::mutex> l(cvmutex);
            workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() { return workQueue.try_pop(miniBatch) or doneParsing; });
        }
        if (miniBatch != nullptr) {
        //while(workQueue.try_pop(miniBatch)) {
            ++activeBatches;
            size_t batchReads{0};
            double logForgettingMass = miniBatch->logForgettingMass;
            std::vector<AlignmentGroup<FragT>*>& alignmentGroups = *(miniBatch->alignments);

            using TranscriptID = size_t;
            using HitIDVector = std::vector<size_t>;
            using HitProbVector = std::vector<double>;

            std::unordered_map<TranscriptID, std::vector<FragT*>> hitList;
            for (auto alnGroup : alignmentGroups) {
                for (auto& a : alnGroup->alignments()) {
                    auto transcriptID = a.transcriptID();
                    if (transcriptID < 0 or transcriptID >= refs.size()) {
                        std::cerr << "Invalid Transcript ID: " << transcriptID << "\n";
                    }
                    hitList[transcriptID].emplace_back(&a);
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
                    auto firstTranscriptID = alnGroup->alignments().front().transcriptID();
                    std::unordered_set<size_t> observedTranscripts;
                    for (auto& aln : alnGroup->alignments()) {
                        auto transcriptID = aln.transcriptID();
                        auto& transcript = refs[transcriptID];
                        transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                        double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;

                        double logFragProb = sailfish::math::LOG_1;
                        //double fragLength = aln.fragLen();

                        switch (aln.fragType()) {
                            case ReadType::SINGLE_END:
                                if (aln.isLeft() and transcript.RefLength - aln.left() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(transcript.RefLength - aln.left());
                                } else if (aln.isRight() and aln.right() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(aln.right());
                                }
                                break;
                            case ReadType::PAIRED_END:
                                logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln.fragLen()));
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
                        double qualProb = -logRefLength + logFragProb + aln.logQualProb();
                        double transcriptLogCount = transcript.mass(initialRound);

                        if ( transcriptLogCount != LOG_0 ) {

                            double errLike = sailfish::math::LOG_1;
                            if (burnedIn) {
                                //errLike = errMod.logLikelihood(aln, transcript);
                            }

                            //aln.logProb = (transcriptLogCount - logRefLength) + qualProb + errLike;
                            aln.logProb = transcriptLogCount + qualProb + errLike;

                            sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);
                            if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                                refs[transcriptID].addTotalCount(1);
                                observedTranscripts.insert(transcriptID);
                            }
                        } else {
                            aln.logProb = LOG_0;
                        }
                    }
                    // normalize the hits
                    if (sumOfAlignProbs == LOG_0) { std::cerr << "0 probability fragment; skipping\n"; continue; }
                    for (auto& aln : alnGroup->alignments()) {
                        aln.logProb -= sumOfAlignProbs;
                        double r = uni(eng);
                        if (!burnedIn and r < std::exp(aln.logProb)) {
                            //auto transcriptID = aln.transcriptID();
                            //auto& transcript = refs[transcriptID];
                            //errMod.update(aln, transcript, aln.logProb, logForgettingMass);
                            if (aln.fragType() == ReadType::PAIRED_END) {
                                double fragLength = aln.fragLen();//std::abs(aln.read1->core.pos - aln.read2->core.pos) + aln.read2->core.l_qseq;
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
                        clusterForest.mergeClusters<FragT>(alnGroup->alignments().begin(),
                                                           alnGroup->alignments().end());
                        clusterForest.updateCluster(alnGroup->alignments().front().transcriptID(),
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
            /*
            if (std::abs(individualTotal - clustTotal) > 1e-5) {
                std::cerr << "Cluster total = " << clustTotal << ", individual total = " << individualTotal << "\n";
            }
            */
            miniBatch->release(alignmentStructureQueue, alignmentGroupQueue);
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
void quantifyLibrary(
        AlignmentLibrary<FragT>& alnLib,
        bfs::path& outputFile,
        size_t numRequiredFragments,
        uint32_t numThreads) {

    bool burnedIn{false};
    ErrorModel errMod(1.00);
    auto& refs = alnLib.transcripts();
    size_t numTranscripts = refs.size();
    size_t miniBatchSize{1000};
    size_t numObservedFragments{0};

    ClusterForest clusterForest(numTranscripts, refs);

    MiniBatchQueue<AlignmentGroup<FragT>> workQueue;
    size_t maxFragLen = 800;
    size_t meanFragLen = 200;
    size_t fragLenStd = 80;
    size_t fragLenKernelN = 4;
    double fragLenKernelP = 0.5;
    FragmentLengthDistribution fragLengthDist(1.0, maxFragLen,
                                              meanFragLen, fragLenStd,
                                              fragLenKernelN,
                                              fragLenKernelP, 1);
    double logForgettingMass{std::log(1.0)};
    double forgettingFactor{0.60};
    size_t batchNum{0};
    double numMappedReads{0.0};
    bool initialRound{true};

    while (numObservedFragments < numRequiredFragments) {
        if (!initialRound) {
            if (!alnLib.reset()) {
                fmt::print(stderr,
                  "\n\n======== WARNING ========\n"
                  "The provided alignment file: [{}] "
                  "is not a regular file and therefore can't be read from "
                  "more than once.\n\n"
                  "We observed only {} fragments when we wanted at least {}.\n\n"
                  "Please consider re-running Salmon with these alignments "
                  "as a regular file!\n"
                  "==========================\n\n",
                  alnLib.alignmentFile(), numObservedFragments, numRequiredFragments);
                break;
            }
        }
        std::atomic<bool> doneParsing{false};
        std::condition_variable workAvailable;
        std::mutex cvmutex;
        std::vector<std::thread> workers;
        std::atomic<size_t> activeBatches{0};
        std::atomic<size_t> processedReads{0};
        for (uint32_t i = 0; i < numThreads; ++i) {
            workers.emplace_back(processMiniBatch<FragT>,
                    std::ref(workQueue),
                    std::ref(workAvailable), std::ref(cvmutex),
                    std::ref(refs), std::ref(clusterForest),
                    std::ref(doneParsing), std::ref(activeBatches),
                    std::ref(alnLib.alignmentStructureQueue()),
                    std::ref(alnLib.alignmentGroupQueue()),
                    std::ref(errMod),
                    std::ref(fragLengthDist),
                    std::ref(burnedIn),
                    initialRound,
                    std::ref(processedReads));
        }

        size_t numProc{0};

        std::vector<AlignmentGroup<FragT>*>* alignments = new std::vector<AlignmentGroup<FragT>*>;
        alignments->reserve(miniBatchSize);
        AlignmentGroup<FragT>* ag;
        while (alnLib.getAlignmentGroup(ag)) {
            alignments->push_back(ag);
            if (alignments->size() >= miniBatchSize) {
                ++batchNum;
                if (batchNum > 1) {
                    logForgettingMass += forgettingFactor * std::log(static_cast<double>(batchNum-1)) -
                        std::log(std::pow(static_cast<double>(batchNum), forgettingFactor) - 1);
                }

                MiniBatchInfo<AlignmentGroup<FragT>>* mbi =
                    new MiniBatchInfo<AlignmentGroup<FragT>>(batchNum, alignments, logForgettingMass);
                workQueue.push(mbi);
                {
                    std::unique_lock<std::mutex> l(cvmutex);
                    workAvailable.notify_one();
                }
                alignments = new std::vector<AlignmentGroup<FragT>*>;
                alignments->reserve(miniBatchSize);
            }
            if (numProc % 1000000 == 0) {
                const char RESET_COLOR[] = "\x1b[0m";
                char green[] = "\x1b[30m";
                green[3] = '0' + static_cast<char>(fmt::GREEN);
                char red[] = "\x1b[30m";
                red[3] = '0' + static_cast<char>(fmt::RED);
                fmt::print(stderr, "\r\r{}processed{} {} {}reads{}\n", green, red, numProc, green, RESET_COLOR);
                fmt::print(stderr, "log(forgettingMass) = {}\x1b[A", logForgettingMass);
            }
            ++numProc;
        }
        std::cerr << "\n";
        for (auto& aln : *alignments) {
            aln->alignments().clear();
            delete aln;
            aln = nullptr;
        }
        delete alignments;

        doneParsing = true;
        size_t tnum{0};
        for (auto& t : workers) {
            fmt::print(stderr, "killing thread {} . . . ", tnum++);
            {
                std::unique_lock<std::mutex> l(cvmutex);
                workAvailable.notify_all();
            }
            t.join();
            fmt::print(stderr, "done\r\r");
        }
        fmt::print(stderr, "\n");

        initialRound = false;
        numObservedFragments += alnLib.numMappedReads();
        numMappedReads = alnLib.numMappedReads();
        fmt::print(stderr, "# observed = {} / # required = {}\n",
                   numObservedFragments, numRequiredFragments);
    }

    std::cerr << "writing output \n";

    std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(outputFile.c_str(), "w"), std::fclose);

    fmt::print(output.get(), "# Salmon-Alignment v 0.01\n");
    fmt::print(output.get(), "# Name\tLength\tTPM\tFPKM\tNumReads\n");

    const double logBillion = std::log(1000000000.0);
    const double million = 1000000.0;
    const double logNumFragments = std::log(static_cast<double>(numMappedReads));
    auto clusters = clusterForest.getClusters();
    size_t clusterID = 0;
    for(auto cptr : clusters) {
        /*
           std::vector<double> readCounts;
           std::vector<double> posterior;
           */

        double logClusterMass = cptr->logMass();
        double logClusterCount = std::log(static_cast<double>(cptr->numHits()));

        if (logClusterMass == LOG_0) {
            std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
        }

        bool requiresProjection{false};

        auto& members = cptr->members();
        size_t clusterSize{0};
        for (auto transcriptID : members) {
            Transcript& t = refs[transcriptID];
            t.uniqueCounts = t.uniqueCount();
            t.totalCounts = t.totalCount();
            //clusterCount += t.totalCounts;
        }

        for (auto transcriptID : members) {
            Transcript& t = refs[transcriptID];
            double logTranscriptMass = t.mass(false);
            if (logTranscriptMass == LOG_0) {
                t.projectedCounts = 0;
            } else {
                double logClusterFraction = logTranscriptMass - logClusterMass;
                t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
                requiresProjection |= t.projectedCounts > static_cast<double>(t.totalCounts) or
                    t.projectedCounts < static_cast<double>(t.uniqueCounts);
            }
            ++clusterSize;
        }

        if (clusterSize > 1 and requiresProjection) {
            cptr->projectToPolytope(refs);
        }
        ++clusterID;
    }

    auto& transcripts_ = refs;
    double tfracDenom{0.0};
    for (auto& transcript : transcripts_) {
        tfracDenom += (transcript.projectedCounts / numMappedReads) / transcript.RefLength;
    }

    // Now posterior has the transcript fraction
    for (auto& transcript : transcripts_) {
        double logLength = std::log(transcript.RefLength);
        double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
        double count = transcript.projectedCounts;
        //double countTotal = transcripts_[transcriptID].totalCounts;
        //double countUnique = transcripts_[transcriptID].uniqueCounts;
        double fpkm = count > 0 ? fpkmFactor * count : 0.0;
        double npm = (transcript.projectedCounts / numMappedReads);
        double tfrac = (npm / transcript.RefLength) / tfracDenom;
        double tpm = tfrac * million;

        fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
                transcript.RefName, transcript.RefLength,
                tpm, fpkm, count);
    }

    // Write the inferred fragment length distribution
    bfs::path distFileName = outputFile.parent_path() / "flenDist.txt";
    {
        std::unique_ptr<std::FILE, int (*)(std::FILE *)> distOut(std::fopen(distFileName.c_str(), "w"), std::fclose);
        fmt::print(distOut.get(), "{}\n", fragLengthDist.toString());
    }
}

int salmonAlignmentQuantify(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    bool noBiasCorrect{false};
    uint32_t numThreads{4};
    size_t requiredObservations{50000000};

    po::options_description generic("Sailfish read-quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("libtype,l", po::value<std::string>(), "Format string describing the library type")
    ("alignments,a", po::value<std::string>(), "input alignment (BAM) file")
    ("targets,t", po::value<std::string>(), "FASTA format file containing target transcripts")
    ("threads,p", po::value<uint32_t>(&numThreads)->default_value(6), "The number of threads to use concurrently.\n"
                                            "The alignment-based quantification mode of salmon is usually I/O bound\n"
                                            "so until there is a faster multi-threaded SAM/BAM parser to feed the \n"
                                            "quantification threads, one should not expect much of a speed-up beyond \n"
                                            "~6 threads.")
    ("output,o", po::value<std::string>(), "Output quantification directory")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
    ("num_required_obs,n", po::value(&requiredObservations)->default_value(50000000),
                                        "The minimum number of observations (mapped reads) that must be observed before\n"
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the \n"
                                        "input file, then it will be read through multiple times.")
    ("gene_map,g", po::value<std::string>(), "File containing a mapping of transcripts to genes.  If this file is provided\n"
                                        "Sailfish will output both quant.sf and quant.genes.sf files, where the latter\n"
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping\n"
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format\n"
                                        "where each line contains the name of a transcript and the gene to which it belongs\n"
                                        "separated by a tab.  The extension of the file is used to determine how the file\n"
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF\n"
                                        "format; files with any other extension are assumed to be in the simple format");

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

        for (auto& opt : orderedOptions.options) {
            std::cerr << "[ " << opt.string_key << "] => {";
            for (auto& val : opt.value) {
                std::cerr << " " << val;
            }
            std::cerr << " }\n";
        }

        // Verify the gene_map before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("gene_map")) {
            // Make sure the provided file exists
            geneMapPath = vm["gene_map"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                std::cerr << "Could not fine transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'gene_map\' option or provide a valid file\n";
                std::exit(1);
            }
        }

        bfs::path alignmentFile(vm["alignments"].as<std::string>());
        std::cerr << vm["alignments"].as<std::string>() << "!!!!\n";
        if (!bfs::exists(alignmentFile)) {
            std::stringstream ss;
            ss << "The provided alignment file: " << alignmentFile <<
                " does not exist!\n";
            throw std::invalid_argument(ss.str());
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
        // If we made it this far, the output directory exists
        bfs::path outputFile = outputDirectory / "quant.sf";
        // The transcript file contains the target sequences
        bfs::path transcriptFile(vm["targets"].as<std::string>());

        switch (libFmt.type) {
            case ReadType::SINGLE_END:
                {
                AlignmentLibrary<UnpairedRead> alnLib(alignmentFile, transcriptFile, libFmt);
                quantifyLibrary<UnpairedRead>(alnLib, outputFile, requiredObservations, numThreads);
                }
                break;
            case ReadType::PAIRED_END:
                {
                AlignmentLibrary<ReadPair> alnLib(alignmentFile, transcriptFile, libFmt);
                quantifyLibrary<ReadPair>(alnLib, outputFile, requiredObservations, numThreads);
                }
                break;
            default:
                std::cerr << "Cannot quantify library of unknown " <<
                        "format " << libFmt << "\n";
                std::exit(1);
        }
        /*
        if (!noBiasCorrect) {
            // Estimated read length
            double estimatedReadLength = 36;
            // Number of k-mers per read
            double kmersPerRead = 1.0;
            // Total number of mapped kmers
            uint64_t mappedKmers= totalHits;

            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasFeatPath = indexDirectory / "bias_feats.txt";
            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrection(biasFeatPath, estFilePath, estimatedReadLength, kmersPerRead, mappedKmers,
                    20, biasCorrectedFile, nbThreads);

        }
        */

        // Not ready for read / alignment-based bias correction yet
        if (!noBiasCorrect) {
            fmt::print(stderr, "Post-hoc bias correction is not yet supported in salmon; disabling\n");
            noBiasCorrect = true;
        }

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("gene_map")) {
            try {
                sailfish::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory,
                                                            !noBiasCorrect);
            } catch (std::exception& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                                   "estimates. The gene-level file(s) may not exist",
                                   e.what());
            }
        }

    } catch (po::error& e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "============\n";
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << "============\n";
        std::cerr << argv[0] << " read-quant was invoked improperly.\n";
        std::cerr << "For usage information, " <<
            "try " << argv[0] << " read-quant --help\nExiting.\n";
    }
    return 0;
}
