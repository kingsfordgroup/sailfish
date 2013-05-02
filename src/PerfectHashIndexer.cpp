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

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

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

#include "CountDBNew.hpp"
#include "cmph.h"

#include "PerfectHashIndex.hpp"

void buildPerfectHashIndex(std::vector<uint64_t>& keys, std::vector<uint32_t>& counts, 
                           size_t merLen, const std::string& indexBaseName ) {

    size_t nkeys = keys.size();

    std::vector<uint64_t> orderedMers(nkeys, 0);

    FILE* mphf_fd = fopen("temp.mph", "w");

    // Source of keys
    cmph_io_adapter_t *source = cmph_io_struct_vector_adapter((void *)&keys[0],(cmph_uint32)sizeof(uint64_t), 
        0, sizeof(uint64_t), nkeys);

    std::cerr << "reading from JFHash into perfect hash\n";
    cmph_t *hash = nullptr;
    size_t i = 0;
    { 
      boost::timer::auto_cpu_timer t;     
      //Create minimal perfect hash function using the brz algorithm.
      cmph_config_t *config = cmph_config_new(source);
      cmph_config_set_algo(config, CMPH_BDZ);
      cmph_config_set_mphf_fd(config, mphf_fd);
      cmph_config_set_memory_availability(config, 1000);
      cmph_config_set_b(config, 3);
      cmph_config_set_keys_per_bin(config, 1);

      hash = cmph_new(config);
      cmph_config_destroy(config);
      cmph_dump(hash, mphf_fd); 
      cmph_destroy(hash);   
      fclose(mphf_fd);
    }

    //Find key
    mphf_fd = fopen("temp.mph", "r");
    hash = cmph_load(mphf_fd);
    fclose(mphf_fd);

    std::unique_ptr<cmph_t, std::function<void(cmph_t*)>> ownedHash(hash, cmph_destroy);

    std::cerr << "saving keys in perfect hash . . .";
    auto start = std::chrono::steady_clock::now();
      {
        boost::timer::auto_cpu_timer t;
        tbb::parallel_for_each( keys.begin(), keys.end(), 
            [&ownedHash, &orderedMers]( uint64_t k ) {
                char *key = (char*)(&k);
                unsigned int id = cmph_search(ownedHash.get(), key, sizeof(uint64_t));
                orderedMers[id] = k;
            }
        );
        
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
            }
      );
      std::cerr << "writing transcript counts to file " << indexBaseName+".sfc\n";
      auto dthread2 = std::thread( [&thash, indexBaseName]() -> void { thash.dumpCountsToFile(indexBaseName+".sfc"); } );

      dthread1.join();
      std::cerr << "done writing index\n";
      dthread2.join();
      std::cerr << "done writing transcript counts\n";
}


int mainIndex( int argc, char *argv[] ) {
    using std::string;
    namespace po = boost::program_options;

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("thash,t", po::value<string>(), "transcript hash file [Jellyfish format]")
    ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
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

        string thashFile = vm["thash"].as<string>();

        /**
        *  Read in the Jellyfish hash of the transcripts
        */
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

        tbb::task_scheduler_init init(12);
        auto it = transcriptHash.iterator_all();
        size_t i = 0;
        while ( it.next() ) {
            keys[i] = it.get_key();
            counts[i] = it.get_val();
            ++i;
        }

        string phfile = vm["index"].as<string>();
        buildPerfectHashIndex(keys, counts, merLen, phfile);

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    }

}