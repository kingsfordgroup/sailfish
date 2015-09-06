#include <algorithm>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <thread>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

// TBB include
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

//#include "BiasIndex.hpp"
#include "VersionChecker.hpp"
#include "SailfishConfig.hpp"
#include "SailfishUtils.hpp"
#include "SailfishIndex.hpp"
#include "GenomicFeature.hpp"
#include "TranscriptGeneMap.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "ReadLibrary.hpp"
#include "RapMapUtils.hpp"
#include "HitManager.hpp"
#include "SASearcher.hpp"
#include "SACollector.hpp"
#include "EmpiricalDistribution.hpp"

#include "spdlog/spdlog.h"

/****** QUASI MAPPING DECLARATIONS *********/
using MateStatus = rapmap::utils::MateStatus;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
/****** QUASI MAPPING DECLARATIONS  *******/

/****** Parser aliases ***/
using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
/****** Parser aliases ***/


using FragLengthCountMap = std::unordered_map<uint32_t, uint64_t>;

using std::string;

constexpr uint32_t readGroupSize{1000};

int performBiasCorrectionSalmon(boost::filesystem::path featPath,
                          boost::filesystem::path expPath,
                        boost::filesystem::path outPath,
                          size_t numThreads);



void processReadsQuasi(paired_parser* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               const SailfishOpts& sfOpts,
               FragLengthCountMap& flMap,
	           std::mutex& iomutex) {

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};

  bool tooManyHits{false};
  size_t maxNumHits{sfOpts.maxReadOccs};
  size_t readLen{0};

  auto& numObservedFragments = readExp.numObservedFragmentsAtomic();
  auto& validHits = readExp.numMappedFragmentsAtomic();
  auto&totalHits = readExp.numFragHitsAtomic();
  auto& upperBoundHits = readExp.upperBoundHitsAtomic();
  auto& eqBuilder = readExp.equivalenceClassBuilder();

  auto sidx = readExp.getIndex();
  SACollector hitCollector(sidx->quasiIndex());
  SASearcher saSearcher(sidx->quasiIndex());
  rapmap::utils::HitCounters hctr;

  std::vector<QuasiAlignment> leftHits;
  std::vector<QuasiAlignment> rightHits;
  std::vector<QuasiAlignment> jointHits;

  std::vector<uint32_t> txpIDs;
  std::vector<double> auxProbs;
  size_t txpIDsHash{0};

  while(true) {
    typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
        readLen = j->data[i].first.seq.length();
        tooManyHits = false;
        jointHits.clear();
        leftHits.clear();
        rightHits.clear();
        txpIDs.clear();
        auxProbs.clear();
        txpIDsHash = 0;

        bool lh = hitCollector(j->data[i].first.seq,
                               leftHits, saSearcher,
                               MateStatus::PAIRED_END_LEFT);

        bool rh = hitCollector(j->data[i].second.seq,
                               rightHits, saSearcher,
                               MateStatus::PAIRED_END_RIGHT);

        rapmap::utils::mergeLeftRightHits(
                               leftHits, rightHits, jointHits,
                               readLen, maxNumHits, tooManyHits, hctr);

        upperBoundHits += (jointHits.size() > 0);

        if (jointHits.size() > sfOpts.maxReadOccs ) { jointHits.clear(); }

        if (!sfOpts.allowOrphans) {
            if (jointHits.size() > 0 and jointHits.front().mateStatus != MateStatus::PAIRED_END_PAIRED) {
                jointHits.clear();
            }
        }

        if (jointHits.size() > 0) {
            // Are the jointHits paired-end quasi-mappings or orphans?
            bool isPaired = jointHits.front().mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED;

            // This is a unique hit
            if (jointHits.size() == 1 and isPaired) {
                flMap[jointHits.front().fragLen]++;
            }

            // If these aren't paired-end reads --- so that
            // we have orphans --- make sure we sort the
            // mappings so that they are in transcript order
            if (!isPaired) {
                 // Find the end of the hits for the left read
                 auto leftHitEndIt = std::partition_point(
                        jointHits.begin(), jointHits.end(),
                        [](const QuasiAlignment& q) -> bool {
                        return q.mateStatus == rapmap::utils::MateStatus::PAIRED_END_LEFT;
                        });
                 // Merge the hits so that the entire list is in order
                 // by transcript ID.
                 std::inplace_merge(jointHits.begin(), leftHitEndIt, jointHits.end(),
                         [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                            return a.transcriptID() < b.transcriptID();
                         });
            }

            auto auxProb = 1.0 / jointHits.size();
            for (auto& h : jointHits) {
                auto transcriptID = h.transcriptID();
                txpIDs.push_back(transcriptID);
                auxProbs.push_back(auxProb);
                boost::hash_combine(txpIDsHash, transcriptID);
            }
            TranscriptGroup tg(txpIDs, txpIDsHash);
            eqBuilder.addGroup(std::move(tg), auxProbs);
        }

        validHits += (jointHits.size() > 0);
        totalHits += jointHits.size();
        locRead++;
        ++numObservedFragments;
        if (numObservedFragments % 500000 == 0) {
    	    iomutex.lock();
            fmt::print(stderr, "\033[A\r\rprocessed {} fragments\n", numObservedFragments);
            fmt::print(stderr, "hits per frag:  {}", totalHits / static_cast<float>(prevObservedFrags));
            iomutex.unlock();
        }

    } // end for i < j->nb_filled
    prevObservedFrags = numObservedFragments;
  }
}

// SINGLE END
void processReadsQuasi(single_parser* parser,
        ReadExperiment& readExp,
        ReadLibrary& rl,
        const SailfishOpts& sfOpts,
        std::mutex& iomutex) {

    uint64_t prevObservedFrags{1};

    size_t locRead{0};
    uint64_t localUpperBoundHits{0};

    bool tooManyHits{false};
    size_t readLen{0};
    size_t maxNumHits{sfOpts.maxReadOccs};

    auto& numObservedFragments = readExp.numObservedFragmentsAtomic();
    auto& validHits = readExp.numMappedFragmentsAtomic();
    auto&totalHits = readExp.numFragHitsAtomic();
    auto& upperBoundHits = readExp.upperBoundHitsAtomic();
    auto& eqBuilder = readExp.equivalenceClassBuilder();

    auto sidx = readExp.getIndex();
    SACollector hitCollector(sidx->quasiIndex());
    SASearcher saSearcher(sidx->quasiIndex());
    rapmap::utils::HitCounters hctr;
    std::vector<QuasiAlignment> jointHits;

    std::vector<uint32_t> txpIDs;
    std::vector<double> auxProbs;
    size_t txpIDsHash{0};

    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
            readLen = j->data[i].seq.length();
            tooManyHits = false;
            localUpperBoundHits = 0;
            jointHits.clear();
            txpIDs.clear();
            auxProbs.clear();
            txpIDsHash = 0;

            bool lh = hitCollector(j->data[i].seq,
                    jointHits, saSearcher,
                    MateStatus::SINGLE_END);

            upperBoundHits += (jointHits.size() > 0);

            // If the read mapped to > maxReadOccs places, discard it
            if (jointHits.size() > sfOpts.maxReadOccs ) { jointHits.clear(); }

            if (jointHits.size() > 0) {
                auto auxProb = 1.0 / jointHits.size();
                for (auto& h : jointHits) {
                    auto transcriptID = h.transcriptID();
                    txpIDs.push_back(transcriptID);
                    auxProbs.push_back(auxProb);
                    boost::hash_combine(txpIDsHash, transcriptID);
                }
                TranscriptGroup tg(txpIDs, txpIDsHash);
                eqBuilder.addGroup(std::move(tg), auxProbs);
            }

            validHits += (jointHits.size() > 0);
            totalHits += jointHits.size();
            locRead++;
            ++numObservedFragments;
            if (numObservedFragments % 500000 == 0) {
                iomutex.lock();
                fmt::print(stderr, "\033[A\r\rprocessed {} fragments\n", numObservedFragments);
                fmt::print(stderr, "hits per frag:  {}", totalHits / static_cast<float>(prevObservedFrags));
                iomutex.unlock();
            }

        } // end for i < j->nb_filled

        prevObservedFrags = numObservedFragments;
    }
}

void setEffectiveLengthsTrivial(ReadExperiment& readExp,
        const SailfishOpts& sfOpts) {
        auto& transcripts = readExp.transcripts();
        for(size_t txpID = 0; txpID < transcripts.size(); ++txpID) {
            auto& txp = transcripts[txpID];
            double refLen = txp.RefLength;
            if (txp.RefLength <= sfOpts.fragLenDistPriorMean) {
                txp.EffectiveLength = refLen;
            } else {
                // Maybe convolve this with the normal given the variance
                // provided by the user.
                txp.EffectiveLength = refLen - sfOpts.fragLenDistPriorMean;
            }
        }
}


void quasiMapReads(
        ReadExperiment& readExp,
        const SailfishOpts& sfOpts,
        std::mutex& iomutex){

    std::vector<std::thread> threads;
    auto& rl = readExp.readLibraries().front();
    rl.checkValid();

    auto numThreads = sfOpts.numThreads;

    std::unique_ptr<paired_parser> pairedParserPtr{nullptr};
    std::unique_ptr<single_parser> singleParserPtr{nullptr};

    // Remember the fragment lengths that we see in each thread
    std::vector<FragLengthCountMap> flMaps(numThreads);


    // If the read library is paired-end
    // ------ Paired-end --------
    if (rl.format().type == ReadType::PAIRED_END) {

        if (rl.mates1().size() != rl.mates2().size()) {
            sfOpts.jointLog->error("The number of provided files for "
                    "-1 and -2 must be the same!");
            std::exit(1);
        }

        size_t numFiles = rl.mates1().size() + rl.mates2().size();
        char** pairFileList = new char*[numFiles];
        for (size_t i = 0; i < rl.mates1().size(); ++i) {
            pairFileList[2*i] = const_cast<char*>(rl.mates1()[i].c_str());
            pairFileList[2*i+1] = const_cast<char*>(rl.mates2()[i].c_str());
        }

        size_t maxReadGroup{readGroupSize}; // Number of reads in each "job"
        size_t concurrentFile{2}; // Number of files to read simultaneously
        pairedParserPtr.reset(new
                paired_parser(4 * numThreads, maxReadGroup,
                    concurrentFile, pairFileList, pairFileList+numFiles));

        for(int i = 0; i < numThreads; ++i)  {
            // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
            // change value before the lambda below is evaluated --- crazy!
            auto threadFun = [&,i]() -> void {
                processReadsQuasi(
                        pairedParserPtr.get(),
                        readExp,
                        rl,
                        sfOpts,
                        flMaps[i],
                        iomutex);
            };
            threads.emplace_back(threadFun);
        }
        // Join the threads and collect the results from the count maps
        size_t totalObs{0};
        std::map<uint32_t, uint32_t> jointMap;

        for(int i = 0; i < numThreads; ++i) {
            threads[i].join();
            auto& flMap = flMaps[i];
            for (auto& kv : flMap) {
                jointMap[kv.first] += kv.second;
                totalObs += kv.second;
            }
        }

        sfOpts.jointLog->info("Gathered fragment lengths from all threads");

        /** If we have a sufficient number of observations for the empirical
         *  distribution, then use that --- otherwise use the provided prior
         *  mean fragment length.
         **/
        // Note: if "noEffectiveLengthCorrection" is set, so that these values
        // won't matter anyway, then don't bother computing this "expensive"
        // version.
        const size_t numRequiredFLDObs{50000};
        if (totalObs > numRequiredFLDObs and !sfOpts.noEffectiveLengthCorrection) {
            std::vector<uint32_t> vals;
            std::vector<uint32_t> multiplicities;

            vals.reserve(jointMap.size());
            multiplicities.reserve(jointMap.size());

            for (auto& kv : jointMap) {
                vals.push_back(kv.first);
                multiplicities.push_back(kv.second);
            }

            sfOpts.jointLog->info("Building empirical fragment length distribution");
            EmpiricalDistribution empDist(vals, multiplicities);
            sfOpts.jointLog->info("finished building empirical fragment length distribution");
            using BlockedIndexRange =  tbb::blocked_range<size_t>;

            tbb::task_scheduler_init tbbScheduler(sfOpts.numThreads);
            auto& transcripts = readExp.transcripts();

            sfOpts.jointLog->info("Estimating effective lengths");
            sfOpts.jointLog->info("Emp. dist min = {}, Emp. dist max = {}",
                                  empDist.minValue(), empDist.maxValue());

            tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts.size())),
                [&transcripts, &empDist](const BlockedIndexRange& range) -> void {
                    for (auto txpID : boost::irange(range.begin(), range.end())) {
                        auto& txp = transcripts[txpID];
                       /**
                          *  NOTE: Adopted from "est_effective_length" at
                          *  (https://github.com/adarob/eXpress/blob/master/src/targets.cpp)
                          *  originally written by Adam Roberts.
                          */
                        double refLen = txp.RefLength;
                        if (refLen <= empDist.median()) {
                            txp.EffectiveLength = refLen;
                        } else {
                           uint32_t mval = empDist.maxValue();
                           double effectiveLength = 0.0;
                           for (size_t l = empDist.minValue(); l <= std::min(txp.RefLength, mval); ++l) {
                                effectiveLength += empDist.pdf(l) * (txp.RefLength - l + 1.0);
                            }
                            txp.EffectiveLength = effectiveLength;
                        }
                    }
                });
        } else {
            sfOpts.jointLog->warn("Sailfish saw fewer then {} uniquely mapped reads "
                                  "so {} will be used as the mean fragment length for "
                                  "effective length correction", numRequiredFLDObs,
                                  sfOpts.fragLenDistPriorMean);
            setEffectiveLengthsTrivial(readExp, sfOpts);
        }
    } // ------ Single-end --------
    else if (rl.format().type == ReadType::SINGLE_END) {

        char* readFiles[] = { const_cast<char*>(rl.unmated().front().c_str()) };
        size_t maxReadGroup{readGroupSize}; // Number of files to read simultaneously
        size_t concurrentFile{1}; // Number of reads in each "job"
        stream_manager streams( rl.unmated().begin(),
                rl.unmated().end(), concurrentFile);

        singleParserPtr.reset(new single_parser(4 * numThreads,
                    maxReadGroup,
                    concurrentFile,
                    streams));

        for(int i = 0; i < numThreads; ++i)  {
            // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
            // change value before the lambda below is evaluated --- crazy!
            auto threadFun = [&,i]() -> void {
                processReadsQuasi(
                        singleParserPtr.get(),
                        readExp,
                        rl,
                        sfOpts,
                        iomutex);
            };
            threads.emplace_back(threadFun);
        }
        for(int i = 0; i < numThreads; ++i) { threads[i].join(); }
        setEffectiveLengthsTrivial(readExp, sfOpts);
    } // ------ END Single-end --------
}

int mainQuantify(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool biasCorrect{false};
    SailfishOpts sopt;
    sopt.numThreads = std::thread::hardware_concurrency();

    vector<string> unmatedReadFiles;
    vector<string> mate1ReadFiles;
    vector<string> mate2ReadFiles;

    po::options_description generic("\n"
            "basic options");
    generic.add_options()
        ("version,v", "print version string")
        ("help,h", "produce help message")
        ("index,i", po::value<string>()->required(), "Sailfish index")
        ("libType,l", po::value<std::string>()->required(), "Format string describing the library type")
        ("unmatedReads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
         "List of files containing unmated reads of (e.g. single-end reads)")
        ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
         "File containing the #1 mates")
        ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
         "File containing the #2 mates")
        ("allowOrphans", po::bool_switch(&(sopt.allowOrphans))->default_value(true), "Consider orphaned reads as valid hits when "
         "performing lightweight-alignment.  This option will increase sensitivity (allow more reads to map and "
         "more transcripts to be detected), but may decrease specificity as orphaned alignments are more likely "
         "to be spurious.")
        ("threads,p", po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads), "The number of threads to use concurrently.")
        ("output,o", po::value<std::string>()->required(), "Output quantification file.")
        ("geneMap,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
         "Sailfish will output both quant.sf and quant.genes.sf files, where the latter "
         "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
         "should be provided as either a GTF file, or a in a simple tab-delimited format "
         "where each line contains the name of a transcript and the gene to which it belongs "
         "separated by a tab.  The extension of the file is used to determine how the file "
         "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
         "format; files with any other extension are assumed to be in the simple format.")
       ("biasCorrect", po::value(&biasCorrect)->zero_tokens(), "[Experimental]: Output both bias-corrected and non-bias-corrected "
         "qunatification estimates.");

    po::options_description advanced("\n"
            "advanced options");
    advanced.add_options()
        /*
        ("fldMax" , po::value<size_t>(&(sopt.fragLenDistMax))->default_value(800), "The maximum fragment length to consider when building the empirical "
         "distribution")
         */
        ("fldMean", po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(200),
            "If single end reads are being used for quantification, or there are an insufficient "
            "number of uniquely mapping reads when performing paired-end quantification to estimate "
            "the empirical fragment length distribution, then use this value to calculate effective lengths")
        ("fldSD" , po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(80),
            "The standard deviation used in the fragment length distribution for single-end quantification or "
            "when an empirical distribution cannot be learned.")
        ("maxReadOcc,w", po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(200), "Reads \"mapping\" to more than this many places won't be considered.")
        ("noEffectiveLengthCorrection", po::bool_switch(&(sopt.noEffectiveLengthCorrection))->default_value(false), "Disables "
         "effective length correction when computing the probability that a fragment was generated "
         "from a transcript.  If this flag is passed in, the fragment length distribution is not taken "
         "into account when computing this probability.")
        ("useVBOpt", po::bool_switch(&(sopt.useVBOpt))->default_value(false), "Use the Variational Bayesian EM rather than the "
     			"traditional EM algorithm to estimate transcript abundances.");

    po::options_description all("sailfish quant options");
    all.add(generic).add(advanced);

    po::options_description visible("sailfish quant options");
    visible.add(generic).add(advanced);

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(all).run();

        po::store(orderedOptions, vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
                Quant
                ==========
                Perform quasi-mapping-based estimation of
                transcript abundance from RNA-seq reads
                )";
            std::cout << hstring << std::endl;
            std::cout << visible << std::endl;
            std::exit(1);
        }

        po::notify(vm);

        std::stringstream commentStream;
        commentStream << "# sailfish (quasi-mapping-based) v" << sailfish::version << "\n";
        commentStream << "# [ program ] => sailfish \n";
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
                std::cerr << "Could not find transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'geneMap\' option or provide a valid file\n";
                std::exit(1);
            }
        }

        bfs::path outputDirectory(vm["output"].as<std::string>());
        bfs::create_directory(outputDirectory);
        if (!(bfs::exists(outputDirectory) and bfs::is_directory(outputDirectory))) {
            std::cerr << "Couldn't create output directory " << outputDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }

        bfs::path indexDirectory(vm["index"].as<string>());
        bfs::path logDirectory = outputDirectory / "logs";

        sopt.indexDirectory = indexDirectory;
        sopt.outputDirectory = outputDirectory;

        // Create the logger and the logging directory
        bfs::create_directory(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "Logs will be written to " << logDirectory.string() << "\n";

        bfs::path logPath = logDirectory / "sailfish_quant.log";
        // must be a power-of-two
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        std::ofstream logFile(logPath.string());
        if (!logFile.good()) {
            std::cerr << "[WARNING]: Could not open log file --- this seems suspicious!\n";
        }

        auto fileSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(logFile);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("stderrLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        sopt.jointLog = jointLog;
        sopt.fileLog = fileLog;

        // Verify that no inconsistent options were provided
        // {
        // }

        jointLog->info() << "parsing read library format";

        vector<ReadLibrary> readLibraries = sailfish::utils::extractReadLibraries(orderedOptions);

        SailfishIndexVersionInfo versionInfo;
        boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
        versionInfo.load(versionPath);

        ReadExperiment experiment(readLibraries, indexDirectory, sopt);

        // Parameter validation
        // If we're allowing orphans, make sure that the read libraries are paired-end.
        // Otherwise, this option makes no sense.
        if (sopt.allowOrphans) {
            for (auto& rl : readLibraries) {
                if (!rl.isPairedEnd()) {
                    jointLog->error("You cannot specify the --allowOrphans argument "
                            "for single-end libraries; exiting!");
                    std::exit(1);
                }
            }
        }
        // end parameter validation


        // This will be the class in charge of maintaining our
        // rich equivalence classes
        experiment.equivalenceClassBuilder().start();

        std::mutex ioMutex;
        fmt::print(stderr, "\n\n");
        quasiMapReads(experiment, sopt, ioMutex);
        fmt::print(stderr, "Done Quasi-Mapping \n\n");

        experiment.equivalenceClassBuilder().finish();
        // Now that the streaming pass is complete, we have
        // our initial estimates, and our rich equivalence
        // classes.  Perform further optimization until
        // convergence.
        CollapsedEMOptimizer optimizer;
        jointLog->info("Starting optimizer");
        //sailfish::utils::normalizeAlphas(sopt, experiment);
        optimizer.optimize(experiment, sopt, 0.01, 10000);
        jointLog->info("Finished optimizer");

        size_t tnum{0};

        jointLog->info("writing output \n");

        bfs::path estFilePath = outputDirectory / "quant.sf";

        commentStream << "# [ mapping rate ] => { " << experiment.mappingRate() * 100.0 << "\% }\n";
        commentString = commentStream.str();

        sailfish::utils::writeAbundancesFromCollapsed(
                sopt, experiment, estFilePath, commentString);

        /*
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
        */

        if (biasCorrect) {
            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasFeatPath = indexDirectory / "bias_feats.txt";
            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrectionSalmon(biasFeatPath, estFilePath, biasCorrectedFile, sopt.numThreads);
        }

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("geneMap")) {
            try {
                sailfish::utils::generateGeneLevelEstimates(geneMapPath,
                        outputDirectory,
                        biasCorrect);
            } catch (std::invalid_argument& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                        "estimates. The gene-level file(s) may not exist",
                        e.what());
            }
        }

        jointLog->flush();
        logFile.close();

    } catch (po::error &e) {
        std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
        std::exit(1);
    }


    return 0;
}
