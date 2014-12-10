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


#include <boost/thread/thread.hpp>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <cassert>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "jellyfish/config.h"
#include "jellyfish/err.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/jellyfish.hpp"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include "cmph.h"
#include "CountDBNew.hpp"
// #include "LookUpTableUtils.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "PerfectHashIndex.hpp"
#include "spdlog/spdlog.h"

void buildPerfectHashIndex(bool canonical, std::vector<uint64_t>& keys, std::vector<uint32_t>& counts,
                           size_t merLen, const boost::filesystem::path& indexBasePath) {

    namespace bfs = boost::filesystem;
    size_t nkeys = keys.size();

    std::vector<uint64_t> orderedMers(nkeys, 0);

    // Source of keys -- oh C, how I love thee
    cmph_io_adapter_t *source = cmph_io_struct_vector_adapter(static_cast<void *>(&keys[0]),
                                                              static_cast<cmph_uint32>(sizeof(uint64_t)),
                                                              0, sizeof(uint64_t), nkeys);

    std::cerr << "Building a perfect hash with " << nkeys << " keys from the Jellyfish hash.\n";
    cmph_t *hash = nullptr;
    size_t i = 0;
    {
      boost::timer::auto_cpu_timer t;
      //Create minimal perfect hash function using the brz algorithm.
      cmph_config_t *config = cmph_config_new(source);
      cmph_config_set_algo(config, CMPH_BDZ);
      cmph_config_set_memory_availability(config, 1000);
      cmph_config_set_b(config, 3);
      cmph_config_set_keys_per_bin(config, 1);

      hash = cmph_new(config);
    }

    assert(hash != nullptr);
    std::unique_ptr<cmph_t, std::function<void(cmph_t*)>> ownedHash(hash, cmph_destroy);

    std::cerr << "saving keys in perfect hash . . .";
    auto start = std::chrono::steady_clock::now();
    {
      boost::timer::auto_cpu_timer t;
      tbb::parallel_for_each( keys.begin(), keys.end(),
        [&ownedHash, &orderedMers]( uint64_t k ) -> void {
          char *key = (char*)(&k);
          unsigned int id = cmph_search(ownedHash.get(), key, sizeof(uint64_t));
          orderedMers[id] = k;
        });

    }

    std::cerr << "done\n";
    auto end = std::chrono::steady_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cerr << "took: " << static_cast<double>(ms.count()) / keys.size() << " us / key\n";

    PerfectHashIndex phi(orderedMers, ownedHash, merLen, canonical);

    bfs::path transcriptomeIndexPath(indexBasePath); transcriptomeIndexPath /= "transcriptome.sfi";
    std::cerr << "writing index to file " << transcriptomeIndexPath << "\n";
    auto dthread1 = std::thread( [&phi, transcriptomeIndexPath]() -> void { phi.dumpToFile(transcriptomeIndexPath.string()); } );

    auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
    auto phiPtr = std::shared_ptr<PerfectHashIndex>(&phi, del);
    CountDBNew thash( phiPtr );

    tbb::parallel_for( size_t{0}, keys.size(),
      [&thash, &keys, &counts]( size_t idx ) {
        auto k = keys[idx]; auto c = counts[idx];
        thash.inc(k, c);
      });

    bfs::path transcriptomeCountPath(indexBasePath); transcriptomeCountPath /= "transcriptome.sfc";

    std::cerr << "writing transcript counts to file " << transcriptomeCountPath << "\n";
    auto dthread2 = std::thread( [&thash, transcriptomeCountPath]() -> void {
                                  thash.dumpCountsToFile(transcriptomeCountPath.string());
                                });

    dthread1.join();
    std::cerr << "done writing index\n";
    dthread2.join();
    std::cerr << "done writing transcript counts\n";
}

//int count_main(int argc, char* argv[]);
int jellyfish_count_main(int argc, char *argv[]);

int runJellyfish(bool canonical,
                 uint32_t merLen,
                 uint32_t numThreads,
                 const std::string& outputStem,
                 std::vector<std::string>& inputFiles) {

    namespace bfs = boost::filesystem;
    bfs::path jfTimePath(outputStem); jfTimePath /= "jf.time";
    bfs::path jfCountPath(outputStem); jfCountPath /= "jf.counts";
    std::stringstream argStream;
    std::stringstream inputFilesStream;
    if (canonical) { argStream << "-C "; }
    argStream << "--timing=" << jfTimePath.string() << " ";
    argStream << "-m " << merLen << " ";
    argStream << "-t " << numThreads << " ";
    argStream << "-o " <<  jfCountPath.string() << " ";

    uintmax_t sumFileSizes{0};
    for (auto& fn : inputFiles) {
        sumFileSizes += boost::filesystem::file_size(fn);
        std::cerr << "file: " << fn << " has size " << sumFileSizes << " bytes\n";
        inputFilesStream << fn << " ";
    }

    // Set the hash size based on the input file sizes
    argStream << "-s " << std::llround(1.2 * sumFileSizes) << " ";
    argStream << inputFilesStream.str();

    std::string argString = argStream.str();
    argString.pop_back();
    std::cerr << "argString = [" << argString << "]\n";
    //boost::trim(argString);

    std::cerr << "running jellyfish with " << argString << "\n";

    // Run Jellyfish as an separate process.
    // This will force the mmapped memory to be cleaned up.
    auto pid = fork();
    if (pid == 0) { //child
        // Create the argv array for the main call to Jellyfish
        std::vector<std::string> argStrings;
        boost::split(argStrings, argString, boost::is_any_of(" "));

        char ** jfargs = new char*[argStrings.size()];
        for (size_t i : boost::irange({0}, argStrings.size())) {
            jfargs[i] = new char[argStrings[i].size()+1];
            std::strcpy(jfargs[i], argStrings[i].c_str());
        }

        std::cerr << "In Jellyfish process. Counting transcript kmers\n";
        int jfRet = jellyfish_count_main(argStrings.size(), jfargs);
        std::cerr << "Finished call to jellyfish process\n";
        for (size_t i : boost::irange({0}, argStrings.size())) {
            delete jfargs[i];
        }
        delete [] jfargs;
        std::cerr << "Exiting jellyfish thread\n";
        std::exit(jfRet);

    } else if (pid < 0) { // fork failed!

        std::cerr << "FATAL ERROR: Failed to spawn Jellyfish process. Exiting\n";
        std::exit(-1);
    } else { // parent

        int status = -1;
        std::cerr << "waiting on " << pid << "\n";
        auto waitRet = waitpid(pid, &status, 0); // wait on the Jellyfish process
        std::cerr << "waited on " << pid << ", waitpid return was " << waitRet << "\n";
        bool statusOK = (status == 0);
        std::cerr << "Jellyfish status " << ((statusOK) ? "OK" : "NOT OK")
                  << " -- return code " << status << "\n";
        assert(statusOK);
    }
    return 0;
    /*
    // Run Jellyfish as an external process using the shell.
    // This will force the mmapped memory to be cleaned up.
    auto jfFile = popen(argString.c_str(), "r");
    int jfRet = pclose(jfFile);
    std::cerr << "Jellyfish finished with status: " << jfRet << "\n";
    */
}

void buildLUTs(
  const std::vector<std::string>& transcriptFiles, //!< File from which transcripts are read
  PerfectHashIndex& transcriptIndex,               //!< Index of transcript kmers
  CountDBNew& transcriptHash,                      //!< Count of kmers in transcripts
  TranscriptGeneMap& tgmap,                        //!< Transcript => Gene map
  const std::string& tlutfname,                    //!< Transcript lookup table filename
  const std::string& klutfname,                    //!< Kmer lookup table filename
  uint32_t numThreads                              //!< Number of threads to use in parallel
  );

int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    boost::filesystem::path outFilePath,
    bool useStreamingParser,
    size_t numThreads);

int mainIndex( int argc, char *argv[] ) {
    using std::string;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();
    bool useStreamingParser = true;

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<std::vector<string>>()->multitoken()->required(), "Transcript fasta file(s)." )
    ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("kmerSize,k", po::value<uint32_t>()->required(), "Kmer size.")
    ("out,o", po::value<string>()->required(), "Output stem [all files needed by Sailfish will be of the form stem.*].")

      //("canonical,c", po::bool_switch(), "Passing this flag in forces all processing to be done on canonical kmers.\n"
      //                                       "This means transcripts will be mapped to their canonical kmer multiset and\n"
      //                                       "that directionality will not be considered when mapping kmers from reads.\n"
      //                                       "This slightly increases ambiguity in the isoform estimation step, but\n"
      //                                       "is generally faster than non-canonical processing (about twice as fast\n"
      //                                       "when counting kmers in the reads).\n")

    //("thash,t", po::value<string>(), "transcript hash file [Jellyfish format]")
    //("index,i", po::value<string>(), "transcript index file [Sailfish format]")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use concurrently.")
    ("force,f", po::bool_switch(), "" )
    ;

    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
index
==========
Builds a perfect hash-based Sailfish index [index] from
the Jellyfish database [thash] of the transcripts.
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        uint32_t merLen = vm["kmerSize"].as<uint32_t>();
        string outputStem = vm["out"].as<string>();
        std::vector<string> transcriptFiles = vm["transcripts"].as<std::vector<string>>();
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        bool force = vm["force"].as<bool>();
        // temporarily deprecated
        // bool canonical = vm["canonical"].as<bool>();
        bool canonical = false;

        // Check to make sure that the specified output directory either doesn't exist, or is
        // a valid path (e.g. not a file)
        namespace bfs = boost::filesystem;
        if (bfs::exists(outputStem) and !bfs::is_directory(outputStem)) {
            std::cerr << "The provided output path [" << outputStem << "] " <<
                         "already exists and is not a directory\n.";
            std::cerr << "Please either provide a different path or " <<
                         "delete the existing file.\n";
            std::exit(1);
        }

        bool mustRecompute = false;
        if (!bfs::exists(outputStem)) {
            mustRecompute = true;
            try {
                bool success = bfs::create_directory(outputStem);
                if (!success) { throw std::runtime_error("unspecified error creating file."); }
            } catch ( std::exception& e ) {
                std::cerr << "Creation of " << outputStem << " failed [" << e.what() << "]\n.";
                std::cerr << "Exiting.\n";
                std::exit(1);
            }
        }

        bfs::path outputPath(outputStem);

        // create the directory for log files
        bfs::path logDir = outputPath / "logs";
        boost::filesystem::create_directory(logDir);

        bfs::path logPath = logDir / "sailfish_index.log";
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("consoleLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        std::cerr << "writing log to " << logPath.string() << "\n";

        // First, compute the transcript features in case the user
        // ever wants to bias-correct his / her results
        bfs::path transcriptBiasFile(outputPath); transcriptBiasFile /= "bias_feats.txt";

        std::cerr << "computeBiasFeatures( {";
        for (auto& tf : transcriptFiles) {
            std::cerr << "[" << tf << "] ";
        }
        std::cerr << ", " << transcriptBiasFile << ", " << useStreamingParser << ", " << numThreads << ")\n";
        computeBiasFeatures(transcriptFiles, transcriptBiasFile, useStreamingParser, numThreads);

        bfs::path jfHashFile(outputPath); jfHashFile /= "jf.counts";

        mustRecompute = (force or !boost::filesystem::exists(jfHashFile));

        if (!mustRecompute) {
            // Check that the jellyfish has at the given location
            // was computed with the correct kmer length.
            std::cout << "Checking that jellyfish hash is up to date" << std::endl;
        }

        if (mustRecompute) {
            std::cerr << "Running Jellyfish on transcripts\n";
            runJellyfish(canonical, merLen, numThreads, outputStem, transcriptFiles);
            std::cerr << "Jellyfish finished\n";

            bfs::path thashFile = jfHashFile;//vm["thash"].as<string>();
            bool recomputePerfectIndex = mustRecompute;

            // Read in the Jellyfish hash of the transcripts
            std::ifstream transcriptDB(thashFile.c_str());
            if (!transcriptDB.good()) {
                jointLog->error() << "Couldn't open the Jellyfish hash [" << thashFile.string() << "] quitting\n";
                std::exit(-1);
            }
            jellyfish::file_header header;
            header.read(transcriptDB);

            std::cerr << "transcript hash size is " << header.size() << "\n";
            size_t nkeys = header.size();
            // Since JF2, key_len() is in terms of bits so merLen = keyLen / 2;
            size_t merLen = header.key_len() / 2;
            jellyfish::mer_dna::k(merLen);

            std::cerr << "header.key_len = " << merLen << "\n";
            std::vector<uint64_t> keys(nkeys,0);
            std::vector<uint32_t> counts(nkeys,0);

            tbb::task_scheduler_init init(numThreads);
            size_t i = 0;

            if (!header.format().compare(binary_dumper::format)) {
                binary_reader reader(transcriptDB, &header);
                while ( reader.next() ) {
                    keys[i] = reader.key().get_bits(0, 2*merLen);
                    counts[i] = reader.val();
                    ++i;
                }

            } else if (!header.format().compare(text_dumper::format)) {
                text_reader reader(transcriptDB, &header);
                while ( reader.next() ) {
                    keys[i] = reader.key().get_bits(0, 2*merLen);
                    counts[i] = reader.val();
                    ++i;
                }
            } else {
                jointLog->error() << "Unknown Jellyfish hash format. quitting\n";
                std::exit(-1);
            }

            keys.resize(i);
            counts.resize(i);

            bfs::path sfIndexBase(outputPath);
            bfs::path sfIndexFile(sfIndexBase); sfIndexFile /= "transcriptome.sfi";
            buildPerfectHashIndex(canonical, keys, counts, merLen, sfIndexBase);

            TranscriptGeneMap tgmap;
            if (vm.count("tgmap") ) { // if we have a GTF file
                string transcriptGeneMap = vm["tgmap"].as<string>();
                std::cerr << "building transcript to gene map using gtf file [" <<
                             transcriptGeneMap << "] . . .\n";
                auto features = GTFParser::readGTFFile<TranscriptGeneID>(transcriptGeneMap);
                tgmap = sailfish::utils::transcriptToGeneMapFromFeatures( features );
                std::cerr << "done\n";
            } else {
                std::cerr << "building transcript to gene map using transcript fasta file [" <<
                             transcriptFiles[0] << "] . . .\n";
                tgmap = sailfish::utils::transcriptToGeneMapFromFasta(transcriptFiles[0]);
                std::cerr << "there are " << tgmap.numTranscripts() << " transcripts . . . ";
                std::cerr << "done\n";
            }

            { // save transcript <-> gene map to archive
                bfs::path tgmOutPath(outputPath); tgmOutPath /= "transcriptome.tgm";
                std::cerr << "Saving transcritpt to gene map to [" << tgmOutPath << "]\n";
                std::ofstream ofs(tgmOutPath.string(), std::ios::binary);
                boost::archive::binary_oarchive oa(ofs);
                // write class instance to archive
                oa << tgmap;
            } // archive and stream closed when destructors are called

            bfs::path sfIndexPath(outputPath); sfIndexPath /= "transcriptome.sfi";
            std::cerr << "Reading transcript index from [" << sfIndexPath << "] . . .";
            auto sfIndex = PerfectHashIndex::fromFile( sfIndexPath.string() );
            auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
            auto sfIndexPtr = std::shared_ptr<PerfectHashIndex>( &sfIndex, del );
            std::cerr << "done\n";

            bfs::path sfCountPath(outputPath); sfCountPath /= "transcriptome.sfc";
            std::cerr << "Reading transcript counts from [" << sfCountPath.string() << "] . . .";
            auto sfTranscriptCountIndex = CountDBNew::fromFile(sfCountPath.string(), sfIndexPtr);
            std::cerr << "done\n";

            bfs::path tlutPath(outputPath); tlutPath /= "transcriptome.tlut";
            bfs::path klutPath(outputPath); klutPath /= "transcriptome.klut";

            buildLUTs(transcriptFiles, sfIndex, sfTranscriptCountIndex,
                      tgmap, tlutPath.string(), klutPath.string(), numThreads);

        } else {
            std::cerr << "All index files seem up-to-date.\n";
            std::cerr << "To force Sailfish to rebuild the index, use the --force option.\n";
        }

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " index was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }

    return 0;
}
