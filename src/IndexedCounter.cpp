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

#include "boost/timer/timer.hpp"
#include "boost/range/irange.hpp"

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

int mainCount( int argc, char *argv[] ) {

    using std::string;
    namespace po = boost::program_options;

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
    ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
    ("counts,c", po::value<string>(), "File where Sailfish read count is written")
    ("threads,p", po::value<uint32_t>()->default_value(12), "The number of threads to use when counting kmers")
    ;

    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
count
==========
Counts the kmers in the set of reads [reads] which also occur in
the Sailfish index [index].  The resulting set of counts relies on the
same index, and the counts will be written to the file [counts].
)";
            std::cout << hstring <<"\n";
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        string countsFile = vm["counts"].as<string>();

        string sfIndexBase = vm["index"].as<string>();
        string sfTrascriptIndexFile = sfIndexBase+".sfi";

        auto phi = PerfectHashIndex::fromFile(sfTrascriptIndexFile);
        std::cerr << "index contained " << phi.numKeys() << " kmers\n";

        size_t nkeys = phi.numKeys();
        size_t merLen = phi.kmerLength();

        size_t numActors = vm["threads"].as<uint32_t>();
        std::vector<std::thread> threads;

        auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
        auto phiPtr = std::shared_ptr<PerfectHashIndex>(&phi, del);

        std::atomic<uint64_t> readNum{0};
        std::atomic<uint64_t> processedReads{0};

        std::vector<string> readFiles = vm["reads"].as<std::vector<string>>();
        for( auto rf : readFiles ) {
            std::cerr << "readFile: " << rf << ", ";
        }
        std::cerr << "\n";

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

        CountDBNew rhash( phiPtr );

        // Open up the transcript file for reading
        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);

        {
          std::atomic<size_t> unmappedKmers{0};
          boost::timer::auto_cpu_timer t;
          auto start = std::chrono::steady_clock::now();

          // Start the desired number of threads to parse the reads
          // and build our data structure.
          for (size_t k = 0; k < numActors; ++k) {

            threads.push_back( std::thread( 
                [&parser, &readNum, &rhash, &start, &phi, &unmappedKmers, merLen]() -> void {
                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream = parser.new_thread();

                    typedef uint64_t BinMer;
                    std::vector<BinMer> fwdMers;
                    std::vector<BinMer> revMers;

                    while ( (read = stream.next_read()) ) {
                        ++readNum;
                        if (readNum % 500000 == 0) {
                            auto end = std::chrono::steady_clock::now();
                            auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
                            auto nsec = sec.count();
                            auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                            std::cerr << "processed " << readNum << " reads (" << rate << ") reads/s\r\r";
                        }
                        std::string seq = std::string(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
                        auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
                        auto readLen = std::distance( seq.begin(), newEnd );
                        if ( readLen < merLen ) { continue; }
                        size_t numKmers = readLen - merLen + 1;
                        
                        if ( numKmers > fwdMers.size() ) {
                            fwdMers.resize(numKmers);
                            revMers.resize(numKmers);
                        }

                        //size_t offset = 0;

                        // The number of valid hits using the forward strand
                        // and the reverse-complement strand
                        size_t fCount = 0;
                        size_t rCount = 0;
                        auto INVALID = phi.INVALID;

                        for ( auto offset : boost::irange(size_t{0},numKmers) ){
                            auto mer = seq.substr(offset, merLen);
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary(mer.c_str(), merLen);
                            auto rmer = jellyfish::parse_dna::reverse_complement(binMer, merLen);
                            auto binMerId = phi.index(binMer);
                            auto rMerId = phi.index(rmer);
                            fwdMers[offset] = binMer;
                            revMers[offset] = rmer;
                            fCount += (binMerId != INVALID);
                            rCount += (rMerId != INVALID);
                        }


                        //enum Strand : uint32_t { forward, reverse, canonical };
                        //Strand s = (fCount > rCount) ? forward : reverse;

                        auto& mers = (fCount > rCount) ? fwdMers : revMers;
                        for ( auto offset : boost::irange(size_t{0},numKmers) ){
                            bool inserted = rhash.inc(mers[offset]);
                            if (!inserted) { unmappedKmers++; }
                        }

                        // original version
                        // for ( auto offset : boost::irange(0,numKmers) ){
                        //     auto mer = seq.substr( offset, merLen );
                        //     auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
                        //     auto rmer = jellyfish::parse_dna::reverse_complement( binMer, merLen );
                        //     binMer = (binMer < rmer) ? binMer : rmer;
                        //     rhash.inc(binMer);
                        // }

                    }

                }) 
            );

        }

        // Wait for all of the threads to finish
        for ( auto& thread : threads ){ thread.join(); }
        std::cerr << "\n" << std::endl;
        rhash.dumpCountsToFile(countsFile);

        // Total kmers
        size_t totalCount = 0;
        for (auto i : boost::irange(size_t(0), rhash.kmers().size())) {
            totalCount += rhash.atIndex(i);
        }
        std::cerr << "There were " << totalCount << ", kmers; " << unmappedKmers << " could not be mapped\n";
        std::cerr << "Mapped " << 
                     (totalCount / static_cast<double>(totalCount + unmappedKmers)) * 100.0 << "% of the kmers\n";

        }

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    }
}
