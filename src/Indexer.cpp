#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "boost/timer/timer.hpp"

#include "tbb/parallel_sort.h"
#include "tbb/task_scheduler_init.h"

#include "threadpool.hpp"
#include "jellyfish/parse_dna.hpp"
#include "jellyfish/mapped_file.hpp"
#include "jellyfish/parse_read.hpp"
#include "jellyfish/sequence_parser.hpp"
#include "jellyfish/dna_codes.hpp"
#include "jellyfish/compacted_hash.hpp"
#include "jellyfish/mer_counting.hpp"
#include "jellyfish/misc.hpp"

#include "CountDB.hpp"


int indexMain( int argc, char *argv[] ) {
    using std::string;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("processors,p", po::value<uint32_t>(maxThreads), "Number of threads to use in parallel")
    ("thash,t", po::value<string>(), "transcript hash file [Jellyfish format]")
    ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
    ("idx,i", po::value<string>(), "File where Sailfish transcript index is written")
    ("count,c", po::value<string>(), "File where Sailfish read count is written")
    ;

    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            std::cout << "BuildIndex\n";
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        string thashFile = vm["thash"].as<string>();
        std::vector<string> readFiles = vm["reads"].as<std::vector<string>>();

        string sfIndexFile = vm["idx"].as<string>();
        string countFile = vm["count"].as<string>();

        /**
        *  Read in the Jellyfish has of the transcripts
        */
        mapped_file transcriptDB(thashFile.c_str());
        transcriptDB.random().will_need();
        char typeTrans[8];
        memcpy(typeTrans, transcriptDB.base(), sizeof(typeTrans));

        hash_query_t transcriptHash(thashFile.c_str());
        std::cerr << "transcriptHash size is " << transcriptHash.get_distinct() << "\n";

        size_t nkeys = transcriptHash.get_distinct();
        size_t merLen = transcriptHash.get_mer_len();

        std::vector<uint64_t> mers(nkeys, 0);

        // The CountDB needs a shared_ptr to the kmer vector
        // until we have a better solution, do this.
        auto noDeleter = []( std::vector<uint64_t> *p ) {};
        std::shared_ptr<std::vector<uint64_t>> mersPtr( &mers, noDeleter );

        std::string mer;
        uint64_t bmer;
        size_t i = 0;
        auto it = transcriptHash.iterator_all();
        while ( it.next() ) {
            if ( i % 1000000 == 0 ) {
                std::cerr << "i = " << i << "\n";
            }
            mers[i] = it.get_key();
            ++i;
        }


        {
            boost::timer::auto_cpu_timer t;
            std::cerr << "sorting ... ";
            tbb::parallel_sort(mers.begin(), mers.end());
            std::cerr << "done\n";
        }

        CountDB transcriptIndex(mersPtr, merLen);

        tbb::parallel_for( size_t(0), mers.size(), size_t(1),
          [&mers, &transcriptHash, &transcriptIndex]( size_t i ) -> void {
            uint64_t kmer = mers[i];
            transcriptIndex.inc( kmer, transcriptHash[kmer] );
          }
        );

        transcriptIndex.dumpToFile( sfIndexFile );

        CountDB readHash(mersPtr, merLen);

        char** fnames = new char*[readFiles.size()];
        size_t z{0};
        size_t numFnames{0};
        for ( auto& s : readFiles ){
            // Ugly, yes?  But this is not as ugly as the alternatives.
            // The char*'s contained in fname are owned by the readFiles
            // vector and need not be manually freed.
            fnames[numFnames] = const_cast<char*>(s.c_str());
            ++numFnames;
        }

        // Open up the transcript file for reading
        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);
        
        size_t numActors = vm["processors"].as<uint32_t>();
        boost::threadpool::pool tp(numActors);
        //std::vector< std::atomic<uint32_t> > readHash( mers.size() );
        std::atomic<uint64_t> readNum {0};
        {
            boost::timer::auto_cpu_timer t;
            for ( size_t a = 0; a < numActors; ++a ) {
                auto worker = [&]() -> void {
                    jellyfish::parse_read::read_t *read;
                    jellyfish::parse_read::thread stream = parser.new_thread();
                    while ( (read = stream.next_read()) ) {
                        ++readNum;
                        if (readNum % 100000 == 0) {
                            std::cerr << "readNum = " << readNum << "\n";
                        }
                        std::string seq = std::string(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
                        auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
                        auto readLen = std::distance( seq.begin(), newEnd );
                        size_t numKmers = readLen - merLen + 1;
                        size_t offset = 0;
                        while ( offset < numKmers ) {
                            auto mer = seq.substr( offset, merLen );
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
                            auto rmer = jellyfish::parse_dna::reverse_complement( binMer, merLen );
                            binMer = (binMer < rmer) ? binMer : rmer;
                            readHash.inc(binMer);
                            ++offset;
                        }
                    }
                };
                tp.schedule(worker);
            }
            tp.wait();
            std::cerr << "timing for counting all reads\n";
        }

        readHash.dumpCountsToFile( countFile );

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    }

}

