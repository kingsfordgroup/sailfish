/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


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

#if HAVE_LOGGER
#include "g2logworker.h"
#include "g2log.h"
#endif

#include "BiasIndex.hpp"
#include "CountDBNew.hpp"
#include "VersionChecker.hpp"
#include "SailfishConfig.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "TranscriptGeneMap.hpp"
#include "CollapsedIterativeOptimizer.hpp"
#include "LibraryFormat.hpp"
#include "ReadLibrary.hpp"

using std::string;

int performBiasCorrection(boost::filesystem::path featPath,
                          boost::filesystem::path expPath,
                          double estimatedReadLength,
                          double kmersPerRead,
                          uint64_t mappedKmers,
                          uint32_t merLen,
                          boost::filesystem::path outPath,
                          size_t numThreads);

//int mainCount(int argc, char* argv[]);
int mainCount(uint32_t numThreads, const std::string& indexBase,
              /*const LibraryFormat& libFmt,
              const std::vector<string>& undirReadFiles, const std::vector<string>& fwdReadFiles,
              const std::vector<string>& revReadFiles,
              */
              const std::vector<ReadLibrary>& readLibraries,
              const std::string& countFileOut,
              bool discardPolyA);


int runIterativeOptimizer(int argc, char* argv[] ) {
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  string cmdString = "estimate";

  try{

   bool poisson = false;
   bool noBiasCorrect = false;
   double minAbundance{0.01};
   double maxDelta{std::numeric_limits<double>::infinity()};
   uint32_t maxThreads = std::thread::hardware_concurrency();
   size_t numIter;

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ("cfg,f", po::value< string >(), "config file")
    ;

    po::options_description config("Configuration");
    config.add_options()
      //("genes,g", po::value< std::vector<string> >(), "gene sequences")
      ("min_abundance, m", po::value<double>(&minAbundance)->default_value(0.0),
       "transcripts with abundance (KPKM) lower than this will be reported at 0.")
      ("counts,c", po::value<string>(), "count file")
      ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
      ("bias,b", po::value<string>(), "bias index prefix (without .bin/.dict)")
      //("thash,t", po::value<string>(), "transcript jellyfish hash file")
      ("output,o", po::value<string>(), "output file")
      ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
      ("delta,d", po::value<double>(&maxDelta)->default_value(5e-3), "consider the optimization to have converged if the relative change in \n"
       "the estimated abundance of all transcripts is below this threshold")
      //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
      ("filter,f", po::value<double>()->default_value(0.0), "during iterative optimization, remove transcripts with a mean less than filter")
      ("iterations,n", po::value<size_t>(&numIter)->default_value(1000), "number of iterations to run the optimzation")
      ("lutfile,l", po::value<string>(), "Lookup table prefix")
      ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);

    if ( vm.count("version") ) {
      std::cout << "version : " << sailfish::version <<"\n";
      std::exit(0);
    }

    if ( vm.count("help") ){
      std::cout << "Sailfish\n";
      std::cout << programOptions << std::endl;
      std::exit(0);
    }

    if ( vm.count("cfg") ) {
      std::cerr << "have detected configuration file\n";
      string cfgFile = vm["cfg"].as<string>();
      std::cerr << "cfgFile : [" << cfgFile << "]\n";
      po::store(po::parse_config_file<char>(cfgFile.c_str(), programOptions, true), vm);
    }
    po::notify(vm);

    bool computeBiasCorrection = !noBiasCorrect;
    uint32_t numThreads = vm["threads"].as<uint32_t>();
    tbb::task_scheduler_init init(numThreads);

    string hashFile = vm["counts"].as<string>();
    //std::vector<string> genesFile = vm["genes"].as<std::vector<string>>();
    //string transcriptHashFile = vm["thash"].as<string>();

    string sfIndexBase = vm["index"].as<string>();

    bfs::path sfIndexBasePath(vm["index"].as<string>());

    string sfIndexFile = sfIndexBase+".sfi";
    string sfTrascriptCountFile = sfIndexBase+".sfc";
    bfs::path outputFilePath = bfs::path(vm["output"].as<string>());

    bfs::path logDir = outputFilePath.parent_path() / "logs";

#if HAVE_LOGGER
    std::cerr << "writing logs to " << logDir.string() << "\n";
    g2LogWorker logger(argv[0], logDir.string());
    g2::initializeLogging(&logger);
#endif

    double minMean = vm["filter"].as<double>();
    string lutprefix = vm["lutfile"].as<string>();
    auto tlutfname = lutprefix + ".tlut";
    auto klutfname = lutprefix + ".klut";
    auto kmerEquivClassFname = bfs::path(tlutfname);
    kmerEquivClassFname = kmerEquivClassFname.parent_path();
    kmerEquivClassFname /= "kmerEquivClasses.bin";

    TranscriptGeneMap tgm;
    { // read the serialized transcript <-> gene map from file
      string tgmFile = sfIndexBase+".tgm";
      std::cerr << "Reading the transcript <-> gene map from [" <<
                   tgmFile << "]\n";
#if HAVE_LOGGER
      LOG(INFO) << "Read transcript <=> gene map from [" << tgmFile << "]";
#endif
      std::ifstream ifs(tgmFile, std::ios::binary);
      boost::archive::binary_iarchive ia(ifs);
      ia >> tgm;
      std::cerr << "done\n";
    }

    std::cerr << "Reading transcript index from [" << sfIndexFile << "] . . .";
    auto sfIndex = PerfectHashIndex::fromFile( sfIndexFile );
    auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
    auto sfIndexPtr = std::shared_ptr<PerfectHashIndex>( &sfIndex, del );
    std::cerr << "done\n";

    // the READ hash
    std::cerr << "Reading read counts from [" << hashFile << "] . . .";
    auto hash = CountDBNew::fromFile( hashFile, sfIndexPtr );
    std::cerr << "done\n";
    //const std::vector<string>& geneFiles{genesFile};
    auto merLen = sfIndex.kmerLength();

    BiasIndex bidx = vm.count("bias") ? BiasIndex( vm["bias"].as<string>() ) : BiasIndex();

    std::cerr << "Creating optimizer . . .";
    CollapsedIterativeOptimizer<CountDBNew> solver(hash, tgm, bidx, numThreads);
    // IterativeOptimizer<CountDBNew, CountDBNew> solver( hash, transcriptHash, tgm, bidx );
    std::cerr << "done\n";

    std::cerr << "optimizing using iterative optimization [" << numIter << "] iterations";

    // EM
    bool haveCI{false};
    solver.optimize(klutfname, tlutfname, kmerEquivClassFname.string(), numIter, minMean, maxDelta);

    // VB
    //bool haveCI{true};
    //solver.optimizeVB(klutfname, tlutfname, kmerEquivClassFname.string(), numIter, minMean, maxDelta);

    std::stringstream headerLines;
    headerLines << "# [sailfish version]\t" << sailfish::version << "\n";
    headerLines << "# [kmer length]\t" << sfIndex.kmerLength() << "\n";
    headerLines << "# [using canonical kmers]\t" << (sfIndex.canonical() ? "true" : "false") << "\n";
    headerLines << "# [command]\t";
    for (size_t i : boost::irange(size_t(0), static_cast<size_t>(argc))) { headerLines << argv[i] << " "; }
    headerLines << "\n";

    solver.writeAbundances(outputFilePath, headerLines.str(), minAbundance, haveCI);

    bool applyCoverageFilter{false};

    if (computeBiasCorrection) {

        // Estimated read length
        double estimatedReadLength = hash.averageLength();
        // Number of k-mers per read
        double kmersPerRead = ((estimatedReadLength - hash.kmerLength()) + 1);
        // Total number of mapped kmers
        uint64_t mappedKmers= 0;
        for (auto i : boost::irange(size_t(0), hash.kmers().size())) {
            mappedKmers += hash.atIndex(i);
        }

        auto origExpressionFile = outputFilePath;

        sfIndexBasePath.remove_filename();
        outputFilePath.remove_filename();

        auto biasFeatPath = sfIndexBasePath / "bias_feats.txt";
        //auto expressionFilePath = outputFilePath / "quant.sf";
        auto expressionFilePath = origExpressionFile;
        auto biasCorrectedFile = outputFilePath / "quant_bias_corrected.sf";
        std::cerr << "biasFeatPath = " << biasFeatPath << "\n";
        std::cerr << "expressionFilePath = " << expressionFilePath << "\n";
        std::cerr << "biasCorrectedFile = " << biasCorrectedFile << "\n";
        performBiasCorrection(biasFeatPath, expressionFilePath, estimatedReadLength, kmersPerRead, mappedKmers,
                              hash.kmerLength(), biasCorrectedFile, numThreads);

        if (applyCoverageFilter) {
            auto transcriptKmerMapFile = sfIndexBasePath / "transcripts.map";
            auto filteredOutputFile = outputFilePath / "quant_bias_corrected_filtered.sf";
            solver.applyCoverageFilter_(biasCorrectedFile, transcriptKmerMapFile, filteredOutputFile, minAbundance);
        }

    } else {
        if (applyCoverageFilter) {
            sfIndexBasePath.remove_filename();
            auto outputFile = outputFilePath;
            auto transcriptKmerMapFile = sfIndexBasePath / "transcripts.map";
            outputFilePath.remove_filename();
            auto filteredOutputFile = outputFilePath / "quant_filtered.sf";
            solver.applyCoverageFilter_(outputFile, transcriptKmerMapFile, filteredOutputFile, minAbundance);
        }
    }

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " " << cmdString << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " " << cmdString << " --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}

//int runIterativeOptimizer(int argc, char* argv[]);

int runKmerCounter(const std::string& sfCommand,
                   uint32_t numThreads,
                   const std::string& indexBase,
                   const std::vector<ReadLibrary>& readLibraries,
                   /*                   const LibraryFormat& libFmt,
                   const std::vector<string>& undirReadFiles,
                   const std::vector<string>& fwdReadFiles,
                   const std::vector<string>& revReadFiles,
                   */
                   const std::string& countFileOut,
                   bool discardPolyA) {

    /*
    std::stringstream argStream;

    argStream << sfCommand << " ";

    // if we're ignoring polyA/polyT in reads
    if (discardPolyA) { argStream << "--polya" << " "; }
    argStream << "--index " << indexBase << " ";
    argStream << "--counts " << countFileOut << " ";
    argStream << "--threads " << numThreads << " ";

    std::cerr << "Undirected reads:[";
    // Undirected reads
    if (undirReadFiles.size() > 0) {
        argStream << "--reads=";
        for (auto& rfile : undirReadFiles) {
            std::cerr << rfile << " ";
            argStream << rfile << " ";
        }
    }
    std::cerr << "]\n";

    std::cerr << "Sense reads:[";
    // Sense reads
    if (fwdReadFiles.size() > 0) {
      argStream << "--forward=";
      for (auto& rfile : fwdReadFiles) {
          std::cerr << rfile << " ";
          argStream << rfile << " ";
      }
    }
    std::cerr << "]\n";

    std::cerr << "Anti-sense reads:[";
    // Anti-sense reads
    if (revReadFiles.size() > 0) {
      argStream << "--reverse=";
      for (auto& rfile : revReadFiles) {
          std::cerr << rfile << " ";
          argStream << rfile << " ";
      }
    }
    std::cerr << "]\n";

    std::string argString = argStream.str();
    boost::trim(argString);
    */

    // Run the Sailfish counting phase as an separate process.
    // This will force the mmapped memory to be cleaned up.
    auto pid = fork();
    if (pid == 0) { //child
        /*
        // Create the argv array for the main call to Sailfish
        std::vector<std::string> argStrings;
        boost::split(argStrings, argString, boost::is_any_of(" "));

        char ** args = new char*[argStrings.size()];
        for (size_t i : boost::irange({0}, argStrings.size())) {
            // args[i] = new char[argStrings[i].length()+1];
            // std::copy(argStrings[i].begin(), argStrings[i].end(), args[i]);
            // args[i][argStrings[i].length()] = '\0';
            args[i] = const_cast<char*>(argStrings[i].c_str());
        }

        std::cerr << "In Sailfish Counting thread. Counting read kmers\n";
        int ret = mainCount(argStrings.size(), args);
        delete [] args;
       int ret = mainCount(numThreads, indexBase, libFmt, undirReadFiles,
                            fwdReadFiles, revReadFiles, countFileOut, discardPolyA);

        */
       int ret = mainCount(numThreads, indexBase, readLibraries, countFileOut, discardPolyA);
        std::exit(ret);

    } else if (pid < 0) { // fork failed!

        std::cerr << "FATAL ERROR: Failed to spawn Sailfish process. Exiting\n";
        std::exit(-1);
    } else { // parent

        int status = -1;
        waitpid(pid, &status, 0); // wait on the Jellyfish process
        std::cerr << "Sailfish terminated with return code " << status << "\n";
        assert( status == 0 );
    }

    return 0;
}


int runSailfishEstimation(const std::string& sfCommand,
                          uint32_t numThreads,
                          const boost::filesystem::path& countFile,
                          const boost::filesystem::path& indexBasePath,
                          //const std::string& tgmap,
                          size_t iterations,
                          const boost::filesystem::path& lookupTableBase,
                          const boost::filesystem::path& outFilePath,
                          bool noBiasCorrect,
                          double minAbundance,
                          double maxDelta) {

  using std::vector;
  using std::string;
  using std::cerr;

    std::stringstream argStream;
    argStream << sfCommand << " ";
    // argStream << "estimate ";
    if (noBiasCorrect) {
        argStream << "--no_bias_correct ";
    }
    argStream << "--index " << indexBasePath.string() << " ";
    argStream << "--counts " << countFile.string() << " ";
    argStream << "--threads " << numThreads << " ";
    //argStream << "--tgmap " << tgmap << " ";
    argStream << "--lutfile " << lookupTableBase.string() << " ";
    argStream << "--iterations " << iterations << " ";
    argStream << "--min_abundance " << minAbundance << " ";
    argStream << "--delta " << maxDelta << " ";
    argStream << "--out " << outFilePath.string();

    std::string argString = argStream.str();
    boost::trim(argString);

    int ret = 0;
    auto estThread = std::thread(
        [&]() -> void {
            std::vector<std::string> argStrings;
            boost::split(argStrings, argString, boost::is_any_of(" "));

            char ** args = new char*[argStrings.size()];
            for (size_t i : boost::irange({0}, argStrings.size())) {
                args[i] = const_cast<char*>(argStrings[i].c_str());
            }

            std::cerr << "In Sailfish estimation thread. Estimating transcript abundances\n";
            int ret = runIterativeOptimizer(argStrings.size(), args);
            delete [] args;
        }
    );

    /*
    std::cerr << "cmdline: " << argString << "\n";
    auto sfEstFile = popen(argString.c_str(), "r");
    int estRet = pclose(sfEstFile);
    std::cerr << "Sailfish estimate finished with status: " << estRet << "\n";
    return estRet;
    */
   estThread.join();
   return ret;
}

int mainQuantify( int argc, char *argv[] ) {

    using std::vector;
    using std::string;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    string sfCommand = argv[0];
    uint32_t maxThreads = std::thread::hardware_concurrency();
    bool noBiasCorrect = false;
    double minAbundance{0.01};
    double maxDelta;
    size_t iterations;

    vector<string> undirReadFiles;// = vm["reads"].as<std::vector<string>>();
    vector<string> fwdReadFiles;// = vm["forward"].as<std::vector<string>>();
    vector<string> revReadFiles;// = vm["reverse"].as<std::vector<string>>();

    vector<string> unmatedReadFiles;// = vm["reads"].as<std::vector<string>>();
    vector<string> mate1ReadFiles;// = vm["forward"].as<std::vector<string>>();
    vector<string> mate2ReadFiles;// = vm["reverse"].as<std::vector<string>>();



    // format string: (T|TYPE)=(SE|PE):(O|ORIENTATION)=(>>|<>|><):(S|STRAND)=(AS|SA|S|A|U)
    po::options_description generic("Sailfish quant options");
    generic.add_options()
    ("version,v", "print version string.")
    ("help,h", "produce help message.")
    ("index,i", po::value<string>()->required(), "Sailfish index [output of the \"Sailfish index\" command.")
    ("libtype,l", po::value<vector<string>>()->required(), "Format string describing the library type.")
    //("libtype,l", po::value<string>(), "Format string describing the library type")
    ("unmated_reads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
     "List of files containing unmated reads of (e.g. single-end reads).")
    // ("forward,F", po::value<vector<string>>(&fwdReadFiles)->multitoken(),
    //  "List of files containing reads oriented in the \"sense\" direction")
    // ("reverse,R", po::value<vector<string>>(&revReadFiles)->multitoken(),
    //  "List of files containing reads oriented in the \"anti-sense\" direction")
    ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
        "File containing the #1 mates.")
    ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
        "File containing the #2 mates.")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction.")
    ("min_abundance,m", po::value<double>(&minAbundance)->default_value(0.0),
     "transcripts with an abundance (KPKM) lower than this value will be reported at zero.")
    //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("out,o", po::value<string>()->required(), "Basename of file where estimates are written.")
    ("iterations,n", po::value<size_t>(&iterations)->default_value(1000), "number of iterations to run the optimzation.")
    ("delta,d", po::value<double>(&maxDelta)->default_value(5e-3), "consider the optimization to have converged if the relative change in \n"
                                                           "the estimated abundance of all transcripts is below this threshold.")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers.")
    ("force,f", po::bool_switch(), "Force the counting phase to rerun, even if a count databse exists." )
    ("polya,a", po::bool_switch(), "polyA/polyT k-mers should be discarded.")
    ("gene_map,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided\n"
                                        "Sailfish will output both quant.sf and quant.genes.sf files, where the latter\n"
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping\n"
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format\n"
                                        "where each line contains the name of a transcript and the gene to which it belongs\n"
                                        "separated by a tab.  The extension of the file is used to determine how the file\n"
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF\n"
                                        "format; files with any other extension are assumed to be in the simple format.")
    ;

    po::variables_map vm;

    try {
        auto orderedOptions = po::command_line_parser(argc, argv).options(generic).run();
        po::store(orderedOptions, vm);

        if ( vm.count("help") ){
          std::cout << "Sailfish quant\n";
          std::cout << generic << std::endl;
          std::exit(0);
        }

        po::notify(vm);

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

        vector<ReadLibrary> readLibraries = sailfish::utils::extractReadLibraries(orderedOptions);

        /*
        vector<ReadLibrary> readLibraries;
        for (auto& opt : orderedOptions.options) {
            if (opt.string_key == "libtype") {
                LibraryFormat libFmt = sailfish::utils::parseLibraryFormatString(opt.value[0]);
                if (libFmt.check()) {
                    std::cerr << libFmt << "\n";
                } else {
                    std::stringstream ss;
                    ss << libFmt << " is invalid!";
                    throw std::invalid_argument(ss.str());
                }
                readLibraries.emplace_back(libFmt);
            } else if (opt.string_key == "mates1") {
                readLibraries.back().addMates1(opt.value);
            } else if (opt.string_key == "mates2") {
                readLibraries.back().addMates2(opt.value);
            } else if (opt.string_key == "unmated_reads") {
                readLibraries.back().addUnmated(opt.value);
            }
        }
        */

        for (auto& rl : readLibraries) { rl.checkValid(); }
        /*
        // Collect the read libraries
        for (auto& libFmtStr : libFmtStrs) {
            LibraryFormat libFmt = parseLibraryFormatString(libFmtStr);

            if (libFmt.check()) {
                std::cerr << libFmt << "\n";
            } else {
                std::stringstream ss;
                ss << libFmt << " is invalid!";
                throw std::invalid_argument(ss.str());
            }

            // Ensure that we have only unmated reads with a single end library
            if (libFmt.type == ReadType::SINGLE_END) {
                if (nextUnpaired >= unmatedReadFiles.size()) {
                    string e= "You must provide unmated read files with a single-end library type";
                    throw std::invalid_argument(e);
                }
                string unpairedFilename = unmatedReadFiles[nextUnpaired];
                readLibraries.emplace_back(libFmt, unpairedFilename);
                ++nextUnpaired;
            }

            // or #1 and #2 mates with a paired-end library
            if (libFmt.type == ReadType::PAIRED_END) {
                if (nextPairedEnd >= mate1ReadFiles.size() or nextPairedEnd >= mate2ReadFiles.size()) {
                    string e = "You must provide #1 and #2 mated read files with a paired-end library type";
                    throw std::invalid_argument(e);
                }
                string mateOneFilename = mate1ReadFiles[nextPairedEnd];
                string mateTwoFilename = mate2ReadFiles[nextPairedEnd];
                readLibraries.emplace_back(libFmt, mateOneFilename, mateTwoFilename);
                ++nextPairedEnd;
            }

        }
        */
        bfs::path indexBasePath(vm["index"].as<string>());
        bfs::path outputBasePath(vm["out"].as<string>());
        //string tgmap = vm["tgmap"].as<string>();
        //string tgmap = indexBase+".tgm";
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        bool force = vm["force"].as<bool>();
        bool discardPolyA = vm["polya"].as<bool>();

        /*
        ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
        ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
        ("counts,c", po::value<string>(), "File where Sailfish read count is written")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
        */

        if (bfs::exists(outputBasePath) and !bfs::is_directory(outputBasePath)) {
            std::cerr << "The provided output path [" << outputBasePath << "] " <<
                         "already exists and is not a directory\n.";
            std::cerr << "Please either provide a different path or " <<
                         "delete the existing file.\n";
            std::exit(1);
        }

        bool mustRecount = false;
        if (!bfs::exists(outputBasePath)) {
            mustRecount = true;
            try {
                bool success = bfs::create_directory(outputBasePath);
                if (!success) { throw std::runtime_error("unspecified error creating file."); }
            } catch ( std::exception& e ) {
                std::cerr << "Creation of " << outputBasePath << " failed [" << e.what() << "]\n.";
                std::cerr << "Exiting.\n";
                std::exit(1);
            }
        }

        // create the directory for log files
        bfs::path logDir = outputBasePath / "logs";
        boost::filesystem::create_directory(logDir);

        //string countFileOut = outputBase + ".sfc";
        bfs::path countFilePath(outputBasePath); countFilePath /= "reads.sfc";
        bfs::path indexPath(indexBasePath); indexPath /= "transcriptome";

        mustRecount = (force or !boost::filesystem::exists(countFilePath));
        if (mustRecount) {
            //          runKmerCounter(sfCommand, numThreads, indexPath.string(), undirReadFiles, fwdReadFiles, revReadFiles, countFilePath.string(), discardPolyA);
            runKmerCounter(sfCommand, numThreads, indexPath.string(), readLibraries, countFilePath.string(), discardPolyA);

        }


        /*
        //("genes,g", po::value< std::vector<string> >(), "gene sequences")
        ("counts,c", po::value<string>(), "count file")
        ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
        ("bias,b", po::value<string>(), "bias index prefix (without .bin/.dict)")
        //("thash,t", po::value<string>(), "transcript jellyfish hash file")
        ("output,o", po::value<string>(), "output file")
        ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
        ("filter,f", po::value<double>()->default_value(0.0), "during iterative optimization, remove transcripts with a mean less than filter")
        ("iterations,i", po::value<size_t>(), "number of iterations to run the optimzation")
        ("lutfile,l", po::value<string>(), "Lookup table prefix")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
        */

        bfs::path lutBasePath(indexBasePath); lutBasePath /= "transcriptome";
        bfs::path estFilePath(outputBasePath); estFilePath /= "quant.sf";
        runSailfishEstimation(sfCommand, numThreads, countFilePath, indexPath,
                              iterations, lutBasePath, estFilePath,
                              noBiasCorrect, minAbundance, maxDelta);

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("gene_map")) {
           std::cerr << "Computing gene-level abundance estimates\n";
           bfs::path gtfExtension(".gtf");
           auto extension = geneMapPath.extension();

           TranscriptGeneMap tranGeneMap;
           // parse the map as a GTF file
           if (extension == gtfExtension) {
               // Using the custom GTF Parser
                //auto features = GTFParser::readGTFFile<TranscriptGeneID>(geneMapPath.string());
                //tranGeneMap = sailfish::utils::transcriptToGeneMapFromFeatures(features);

               // Using libgff
                tranGeneMap = sailfish::utils::transcriptGeneMapFromGTF(geneMapPath.string(), "gene_id");
           } else { // parse the map as a simple format files
               std::ifstream tgfile(geneMapPath.string());
                tranGeneMap = sailfish::utils::readTranscriptToGeneMap(tgfile);
                tgfile.close();
           }

           std::cerr << "There were " << tranGeneMap.numTranscripts() << " transcripts mapping to "
                     << tranGeneMap.numGenes() << " genes\n";

           sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, estFilePath);
           if (!noBiasCorrect) {
                bfs::path biasCorrectEstFilePath(estFilePath.parent_path());
                biasCorrectEstFilePath /= "quant_bias_corrected.sf";
                sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, biasCorrectEstFilePath);
           }

        }

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "============\n";
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << "============\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
    }

    return 0;
}
