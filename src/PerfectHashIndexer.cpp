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

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include "jellyfish/parse_dna.hpp"
#include "jellyfish/mapped_file.hpp"
#include "jellyfish/parse_read.hpp"
#include "jellyfish/sequence_parser.hpp"
#include "jellyfish/dna_codes.hpp"
#include "jellyfish/compacted_hash.hpp"
#include "jellyfish/mer_counting.hpp"
#include "jellyfish/misc.hpp"

#include "cmph.h"
#include "CountDBNew.hpp"
// #include "LookUpTableUtils.hpp"
#include "utils.hpp"
#include "genomic_feature.hpp"
#include "PerfectHashIndex.hpp"

void buildPerfectHashIndex(std::vector<uint64_t>& keys, std::vector<uint32_t>& counts, 
                           size_t merLen, const std::string& indexBaseName) {

    size_t nkeys = keys.size();

    std::vector<uint64_t> orderedMers(nkeys, 0);

    // Source of keys -- oh C, how I love thee
    cmph_io_adapter_t *source = cmph_io_struct_vector_adapter(static_cast<void *>(&keys[0]),
                                                              static_cast<cmph_uint32>(sizeof(uint64_t)), 
                                                              0, sizeof(uint64_t), nkeys);

    std::cerr << "Building a perfect hash from the Jellyfish hash.\n";
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

    PerfectHashIndex phi(orderedMers, ownedHash, merLen);
      
    std::cerr << "writing index to file " << indexBaseName+".sfi\n";
    auto dthread1 = std::thread( [&phi, indexBaseName]() -> void { phi.dumpToFile(indexBaseName+".sfi"); } );

    auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
    auto phiPtr = std::shared_ptr<PerfectHashIndex>(&phi, del);
    CountDBNew thash( phiPtr );

    tbb::parallel_for( size_t{0}, keys.size(),
      [&thash, &keys, &counts]( size_t idx ) {
        auto k = keys[idx]; auto c = counts[idx];
        thash.inc(k, c);
      });

    std::cerr << "writing transcript counts to file " << indexBaseName+".sfc\n";
    auto dthread2 = std::thread( [&thash, indexBaseName]() -> void { thash.dumpCountsToFile(indexBaseName+".sfc"); } );

    dthread1.join();
    std::cerr << "done writing index\n";
    dthread2.join();
    std::cerr << "done writing transcript counts\n";
}

int count_main(int argc, char* argv[]);

int runJellyfish(uint32_t merLen, 
                 uint32_t numThreads, 
                 const std::string& outputStem, 
                 std::vector<std::string>& inputFiles) {

    size_t hashSize = 8000000000; // eight billion?

    std::stringstream argStream;
    argStream << "jellyfish ";
    argStream << "count ";
    argStream << "--timing=" << outputStem << ".jftime ";
    argStream << "-m " << merLen << " ";
    argStream << "-s " << hashSize << " ";
    argStream << "-t " << numThreads << " ";
    argStream << "-o " << outputStem << ".jfcounts ";
    for (auto& fn : inputFiles) {
        argStream << fn << " ";
    }

    std::string argString = argStream.str();
    boost::trim(argString);

    // Run Jellyfish as an external process.
    // This will force the mmapped memory to be cleaned up.
    auto jfFile = popen(argString.c_str(), "r");
    int jfRet = pclose(jfFile);
    std::cerr << "Jellyfish finished with status: " << jfRet << "\n";

    /*
    std::vector<std::string> argStrings;
    boost::split(argStrings, argString, boost::is_any_of(" "));

    char** args = new char*[argStrings.size()];

    int argc = argStrings.size();
    for (auto i : boost::irange({0}, argc)) {
        args[i] = const_cast<char*>(argStrings[i].c_str());
        std::cerr << "jfargs[" << i << "] = " << args[i] << "\n";
    }

    auto t = std::thread( [=]() -> int { return count_main(argc, args); });
    t.join();
    delete [] args;
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

int mainIndex( int argc, char *argv[] ) {
    using std::string;

    std::cerr << "running indexer\n";
    for (size_t i = 0; i < argc; ++i) {
        std::cerr << "argv[" << i << "] = " << argv[i] << "\n";
    }

    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<std::vector<string>>()->multitoken(), "Transcript fasta file(s)." )
    ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("kmerSize,k", po::value<uint32_t>()->required(), "Kmer size.")
    ("out,o", po::value<string>(), "Output stem [all files needed by Sailfish will be of the form stem.*].")
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
        string transcriptGeneMap = vm["tgmap"].as<string>();

        string jfHashFile = outputStem + ".jfcounts_0";
        bool mustRecompute = (force or !boost::filesystem::exists(jfHashFile));

        if (!mustRecompute) {
            // Check that the jellyfish has at the given location 
            // was computed with the correct kmer length.
            std::cout << "Checking that jellyfish hash is up to date" << std::endl;
        }

        if (mustRecompute) {
            std::cerr << "Running Jellyfish on transcripts\n";
            runJellyfish(merLen, numThreads, outputStem, transcriptFiles);

            std::cerr << "Jellyfish finished\n";

            string thashFile = jfHashFile;//vm["thash"].as<string>();
            bool recomputePerfectIndex = mustRecompute;

            // Read in the Jellyfish hash of the transcripts
            mapped_file transcriptDB(thashFile.c_str());
            transcriptDB.random().will_need();
            char typeTrans[8];
            memcpy(typeTrans, transcriptDB.base(), sizeof(typeTrans));

            hash_query_t transcriptHash(thashFile.c_str());
            std::cerr << "transcriptHash size is " << transcriptHash.get_distinct() << "\n";
            size_t nkeys = transcriptHash.get_distinct();
            size_t merLen = transcriptHash.get_mer_len();

            std::vector<uint64_t> keys(nkeys,0);
            std::vector<uint32_t> counts(nkeys,0);

            tbb::task_scheduler_init init(numThreads);
            auto it = transcriptHash.iterator_all();
            size_t i = 0;
            while ( it.next() ) {
                keys[i] = it.get_key();
                counts[i] = it.get_val();
                ++i;
            }

            string sfIndexBase = outputStem;
            string sfIndexFile = outputStem + ".sfi";
            buildPerfectHashIndex(keys, counts, merLen, sfIndexBase);

            std::cerr << "parsing gtf file [" << transcriptGeneMap  << "] . . . ";
            auto features = GTFParser::readGTFFile<TranscriptGeneID>(transcriptGeneMap);
            std::cerr << "done\n";

            std::cerr << "building transcript to gene map . . .";
            auto tgmap = utils::transcriptToGeneMapFromFeatures( features );
            std::cerr << "done\n";

            std::cerr << "Reading transcript index from [" << sfIndexFile << "] . . .";
            auto sfIndex = PerfectHashIndex::fromFile( sfIndexFile );
            auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
            auto sfIndexPtr = std::shared_ptr<PerfectHashIndex>( &sfIndex, del );
            std::cerr << "done\n";

            string sfTrascriptCountFile = outputStem + ".sfc";
            std::cerr << "Reading transcript counts from [" << sfTrascriptCountFile << "] . . .";
            auto sfTranscriptCountIndex = CountDBNew::fromFile(sfTrascriptCountFile, sfIndexPtr);
            std::cerr << "done\n";

            string tlutfname = outputStem + ".tlut";
            string klutfname = outputStem + ".klut";

            buildLUTs(transcriptFiles, sfIndex, sfTranscriptCountIndex, tgmap, tlutfname, klutfname, numThreads);

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