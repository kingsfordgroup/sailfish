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

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Sailfish count options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
    ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
    ("counts,c", po::value<string>(), "File where Sailfish read count is written")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
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

        std::cerr << "reading index . . . ";
        auto phi = PerfectHashIndex::fromFile(sfTrascriptIndexFile);
        std::cerr << "done\n";
        std::cerr << "index contained " << phi.numKeys() << " kmers\n";

        size_t nkeys = phi.numKeys();
        size_t merLen = phi.kmerLength();

        size_t numActors = vm["threads"].as<uint32_t>();
        tbb::task_scheduler_init init(numActors);
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
        jellyfish::parse_read parser( fnames, fnames+numFnames, 5000);

        {
          std::atomic<size_t> unmappedKmers{0};
          boost::timer::auto_cpu_timer t(std::cerr);
          auto start = std::chrono::steady_clock::now();

          // Start the desired number of threads to parse the reads
          // and build our data structure.
          for (size_t k = 0; k < numActors; ++k) {

            // threads.push_back( std::thread( 
            //     [&parser, &readNum, &rhash, &start, &phi, &unmappedKmers, merLen]() -> void {
            //         // Each thread gets it's own stream
            //         jellyfish::parse_read::read_t* read;
            //         jellyfish::parse_read::thread stream = parser.new_thread();

            //         typedef uint64_t BinMer;
            //         std::vector<BinMer> fwdMers;
            //         std::vector<BinMer> revMers;

            //         while ( (read = stream.next_read()) ) {
            //             ++readNum;
            //             if (readNum % 500000 == 0) {
            //                 auto end = std::chrono::steady_clock::now();
            //                 auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
            //                 auto nsec = sec.count();
            //                 auto rate = (nsec > 0) ? readNum / sec.count() : 0;
            //                 std::cerr << "processed " << readNum << " reads (" << rate << ") reads/s\r\r";
            //             }
                        
            //             std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );

            //             auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
            //             auto readLen = std::distance( seq.begin(), newEnd );
            //             if ( readLen < merLen ) { continue; }
            //             size_t numKmers = readLen - merLen + 1;
                        
            //             if ( numKmers > fwdMers.size() ) {
            //                 fwdMers.resize(numKmers);
            //                 revMers.resize(numKmers);
            //             }

            //             // The number of valid hits using the forward strand
            //             // and the reverse-complement strand
            //             size_t fCount = 0;
            //             size_t rCount = 0;
            //             auto INVALID = phi.INVALID;

            //             for ( auto offset : boost::irange(size_t{0}, numKmers) ){
            //                 auto mer = seq.substr(offset, merLen);
            //                 auto binMer = jellyfish::parse_dna::mer_string_to_binary(mer.c_str(), merLen);
            //                 auto rmer = jellyfish::parse_dna::reverse_complement(binMer, merLen);
            //                 auto binMerId = phi.index(binMer);
            //                 auto rMerId = phi.index(rmer);
            //                 fwdMers[offset] = binMer;
            //                 revMers[offset] = rmer;
            //                 fCount += (binMerId != INVALID);
            //                 rCount += (rMerId != INVALID);
            //             }

            //             auto& mers = (fCount > rCount) ? fwdMers : revMers;
            //             for ( auto offset : boost::irange(size_t{0}, numKmers) ){
            //                 bool inserted = rhash.inc(mers[offset]);
            //                 if (!inserted) { ++unmappedKmers; }
            //             }

            //             // original version
            //             // for ( auto offset : boost::irange(0,numKmers) ){
            //             //     auto mer = seq.substr( offset, merLen );
            //             //     auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
            //             //     auto rmer = jellyfish::parse_dna::reverse_complement( binMer, merLen );
            //             //     binMer = (binMer < rmer) ? binMer : rmer;
            //             //     rhash.inc(binMer);
            //             // }

            //         }

            //     }) 
            // );
            // 
            // 
            
            /** Guillaume inspired fast parser **/
            threads.push_back( std::thread( 
                [&parser, &readNum, &rhash, &start, &phi, &unmappedKmers, merLen]() -> void {
                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream = parser.new_thread();

                    typedef uint64_t BinMer;
                    std::vector<BinMer> fwdMers;
                    std::vector<BinMer> revMers;

                    BinMer lshift(2 * (merLen - 1));
                    BinMer masq((1UL << (2 * merLen)) - 1);
                    BinMer cmlen, kmer, rkmer;

                    size_t numKmers = 0;
                    size_t offset = 0;
                    size_t fCount = 0; size_t rCount = 0;

                    auto INVALID = phi.INVALID;

                    while ( (read = stream.next_read()) ) {
                        ++readNum;
                        if (readNum % 500000 == 0) {
                            auto end = std::chrono::steady_clock::now();
                            auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
                            auto nsec = sec.count();
                            auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                            std::cerr << "processed " << readNum << " reads (" << rate << ") reads/s\r\r";
                        }

                        // we iterate over the entire read
                        const char         *start = read->seq_s;
                        const char * const  end   = read->seq_e;

                        // reset all of the counts
                        offset = fCount = rCount = numKmers = 0;
                        cmlen = kmer = rkmer = 0;

                        // the maximum number of kmers we'd have to store
                        size_t maxNumKmers = std::distance(start, end);
                        // the read must be at least the kmer length
                        if ( maxNumKmers < merLen ) { continue; }

                        if ( maxNumKmers > fwdMers.size()) {
                            fwdMers.resize(maxNumKmers);
                            revMers.resize(maxNumKmers);
                        }

                        // iterate over the read base-by-base
                        while(start < end) {
                            uint_t     c = jellyfish::dna_codes[static_cast<uint_t>(*start++)];

                            // ***** Potentially consider quality values in the future **** /
                            // const char q = *start++;
                            // if(q < q_thresh)
                            //   c = CODE_RESET;

                            switch(c) {
                                case jellyfish::CODE_IGNORE: break;
                                case jellyfish::CODE_COMMENT:
                                  std::cerr << "ERROR\n";
                                  //report_bad_input(*(start-1));
                                // Fall through
                                case jellyfish::CODE_RESET:
                                  cmlen = kmer = rkmer = 0;
                                  break;

                                default:
                                  // form the new kmer
                                  kmer = ((kmer << 2) & masq) | c;
                                  // the new kmer's reverse complement
                                  rkmer = (rkmer >> 2) | ((0x3 - c) << lshift);
                                  // count if the kmer is valid in the forward and
                                  // reverse directions
                                  if(++cmlen >= merLen) {
                                    cmlen = merLen;
                                    auto binMerId = phi.index(kmer);
                                    auto rMerId = phi.index(rkmer);
                                    fwdMers[offset] = kmer;
                                    revMers[offset] = rkmer;
                                    fCount += (binMerId != INVALID);
                                    rCount += (rMerId != INVALID);
                                    ++offset;
                                    ++numKmers;
                                  } // end if
                            } // end switch
                        } // end read

                        // whichever read direction has more valid kmers is the one we choose
                        auto& mers = (fCount > rCount) ? fwdMers : revMers;
                        // insert the relevant kmers into the count index
                        for ( auto offset : boost::irange(size_t{0}, numKmers) ){
                         bool inserted = rhash.inc(mers[offset]);
                         if (!inserted) { ++unmappedKmers; }
                        }

                } // end parse all reads
            }));

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
        std::cerr << "Program Options Error : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception &e) {
        std::cerr << "ERROR: " << argv[0] << " count invoked improperly.\n";
        std::cerr << "Usage\n";
        std::cerr << "=====\n";
        std::cout << generic << std::endl;
        std::exit(1);
    }

    return 0;
}
