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
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <random>
#include <functional>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/thread/thread.hpp>
#include <boost/filesystem.hpp>

#define _GNU_SOURCE
#include <pthread.h>
#include <sched.h>

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
    namespace bfs = boost::filesystem;

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

        // auto cpuset = new cpu_set_t;
        // CPU_ZERO(cpuset);
        // CPU_SET(0, cpuset);
        // if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), cpuset)) {
        //   std::cerr << "COULD NOT SET PROCESSOR AFFINITY!!!\n";
        //   std::exit(1);
        // }



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

        //phi.will_need(0, numActors+1);
        //rhash.will_need(0, numActors+1);

        // Open up the transcript file for reading
        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 5000);

        {
          std::atomic<size_t> unmappedKmers{0};
          boost::timer::auto_cpu_timer t(std::cerr);
          auto start = std::chrono::steady_clock::now();
          bool canonical = phi.canonical();
          //tbb::concurrent_unordered_set<int> assignedCPUs;

          std::atomic<uint32_t> numPaged{0};

          // Start the desired number of threads to parse the reads
          // and build our data structure.
          for (size_t k = 0; k < numActors; ++k) {
            size_t threadIdx = k;
            /** Guillaume inspired fast parser **/

            // If we're only hashing canonical kmers
            if (canonical) {
                threads.emplace_back(std::thread(
                    [&parser, &readNum, &rhash, &start, &phi, &unmappedKmers, merLen]() -> void {
                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream = parser.new_thread();

                    using BinMer =  uint64_t;

                    BinMer lshift(2 * (merLen - 1));
                    BinMer masq((1UL << (2 * merLen)) - 1);
                    BinMer cmlen, kmer, rkmer;

                    auto INVALID = phi.INVALID;

                    uint64_t localUnmappedKmers{0};

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
                        cmlen = kmer = rkmer = 0;

                        // the maximum number of kmers we'd have to store
                        uint32_t maxNumKmers = std::distance(start, end);

                        // tell the readhash about this read's length
                        rhash.appendLength(maxNumKmers);

                        // the read must be at least the kmer length
                        if ( maxNumKmers < merLen ) { continue; }

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
                                    auto mmer = (kmer < rkmer) ? kmer : rkmer;
                                    if ( phi.index(mmer) != INVALID) { rhash.inc(mmer); } else { ++localUnmappedKmers; }
                                  } // end if
                            } // end switch
                        } // end read
                } // end parse all reads
                unmappedKmers += localUnmappedKmers;
            }));

            } else {
              // If we're hashing kmers in both directions to determine
              // the "direction" of reads.

                enum class MerDirection : std::int8_t { FORWARD = 1, REVERSE = 2, BOTH = 3 };

          //   cpu_set_t* cpuset;
          //   pthread_getaffinity_np(threads.back().native_handle(), sizeof(cpu_set_t), cpuset);
          //   for (size_t cpuID = 0; cpuID < maxThreads; ++cpuID) {
          //     if (CPU_ISSET(cpuID, cpuset)) {
          //       if (CPUMap.find(cpuID) == knownCPUs.end()) {
          //         CPUMap[cpuID] = k;
          //       }
          //     }
          //   }
          // for (auto& CPUThreadPair : CPUMap ) {
          //   phi.will_need(k, CPUMap.size());
          //   rhash.will_need(k, CPUMap.size());
          // }

                threads.emplace_back(std::thread(
                    [&parser, &readNum, &rhash, &start, &phi, &unmappedKmers, &k, &numPaged, threadIdx, merLen, numActors]() -> void {
                      //phi.will_need(threadIdx+1, numActors+1);
                      //rhash.will_need(threadIdx+1, numActors+1);
                      ++numPaged;
                      //while ( numPaged < numActors ) { }
                      //if (threadIdx == numActors - 1) { start = std::chrono::steady_clock::now(); }

                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream{parser.new_thread()};
                    
                    

                    using BinMer = uint64_t;
                    std::vector<BinMer> fwdMers;
                    std::vector<BinMer> revMers;

                    BinMer lshift{2 * (merLen - 1)};
                    BinMer masq{(1UL << (2 * merLen)) - 1};
                    BinMer cmlen, kmer, rkmer;

                    size_t numKmers = 0;
                    size_t numRemaining = 0;
                    size_t fCount = 0; size_t rCount = 0;
                    auto dir = MerDirection::BOTH;

                    auto INVALID = phi.INVALID;

                    uint64_t localUnmappedKmers{0};
                    uint64_t locallyProcessedReads{0};
                    while ( (read = stream.next_read()) ) {
                        ++readNum; ++locallyProcessedReads;
                        if (readNum % 250000 == 0) {
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
                        fCount = rCount = numKmers = 0;
                        cmlen = kmer = rkmer = 0;
                        dir = MerDirection::BOTH;

                        uint32_t readLen = std::distance(start, end);

                        // the maximum number of kmers we'd have to store
                        uint32_t maxNumKmers = (readLen >= merLen) ? readLen - merLen + 1 : 0;
                        numRemaining = maxNumKmers;

                        // tell the readhash about this read's length
                        rhash.appendLength(readLen);

                        // the read must be at least the kmer length
                        if ( maxNumKmers == 0 ) { continue; }

                        if ( maxNumKmers > fwdMers.size()) {
                            fwdMers.resize(maxNumKmers);
                            revMers.resize(maxNumKmers);
                        }

                        size_t binMerId{0};
                        size_t rMerId{0};
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

                                // Fall through
                                case jellyfish::CODE_RESET:
                                  cmlen = kmer = rkmer = 0;
                                  break;

                                default:

                                  kmer = ((kmer << 2) & masq) | c;
                                  rkmer = (rkmer >> 2) | ((0x3 - c) << lshift);
                                  // count if the kmer is valid in the forward and
                                  // reverse directions
                                  if(++cmlen >= merLen) {
                                    cmlen = merLen;

                                    // dispatch on the direction
                                    switch (dir) {
                                       // We're certain that more kmers map in the forward direction
                                       // so we only consider the rest of the read in this direction.
                                       case MerDirection::FORWARD:
                                        // get the index of the forward kmer
                                        binMerId = phi.index(kmer); 
                                        if (binMerId != INVALID) {
                                          rhash.incAtIndex(binMerId);
                                          ++fCount;
                                        }                                         
                                        ++numKmers; --numRemaining;
                                        break;
                                       // end case FORWARD
                                       
                                       // We're certain that more kmers map in the reverse direction
                                       // so we only consider the rest of the read in this direction.
                                       case MerDirection::REVERSE:
                                          // get the index of the forward kmer
                                          rMerId = phi.index(rkmer);
                                          if (rMerId != INVALID) {
                                            rhash.incAtIndex(rMerId);
                                            ++rCount;
                                          }
                                          ++numKmers; --numRemaining;
                                          break;
                                       // end case REVERSE
                                       
                                       case MerDirection::BOTH:
                                          // form the new kmer and it's reverse complement

                                          // Find the index of the forward kmer and determine
                                          // whether or not to count it.
                                          binMerId = phi.index(kmer);
                                          fwdMers[fCount] = binMerId;
                                          fCount += (binMerId != INVALID);

                                          // Find the index of the reverse kmer and determine
                                          // whether or not to count it.
                                          rMerId = phi.index(rkmer);
                                          revMers[rCount] = rMerId;
                                          rCount += (rMerId != INVALID);

                                          ++numKmers; --numRemaining;

                                          // Determine if we need to continue looking at both directions
                                          dir = (fCount > (rCount + numRemaining)) ? MerDirection::FORWARD :
                                                (rCount > (fCount + numRemaining)) ? MerDirection::REVERSE : MerDirection::BOTH;

                                          switch (dir) {
                                            case MerDirection::FORWARD:
                                              for (auto i : boost::irange(size_t(0), fCount)) { rhash.incAtIndex(fwdMers[i]); 
                                              }
                                              break;
                                            case MerDirection::REVERSE:
                                              for (auto i : boost::irange(size_t(0), rCount)) { rhash.incAtIndex(revMers[i]); 
                                              }
                                              break;
                                            default:
                                              break;
                                          }                                        
                                      // end case BOTH
                                      
                                    } // end dirction switch
                                  } // end if
                            } // end switch
                        } // end read

                        uint64_t count{0};
                        switch (dir) {

                          // The same number of things mapped in both directions.  In
                          // this case, we _arbitrarily_ choose the forward kmers. We haven't 
                          // actually incremented counts yet, so we do that here.
                          case MerDirection::BOTH:
                            if (dir == MerDirection::BOTH) {
                              for (auto i : boost::irange(size_t(0), fCount)) { rhash.incAtIndex(fwdMers[i]); 
                              }
                            }
                            count = fCount;
                            break;

                          // More things mapped in the forward direction
                          case MerDirection::FORWARD:
                            count = fCount; break;

                          // More things mapped in the reverse direction                            
                          case MerDirection::REVERSE:
                            count = rCount; break;
                        }

                        // the number of unmapped kmers is just the total kmers in this read
                        // minus the number that mapped.
                        localUnmappedKmers += (numKmers - count);
                        
                } // end parse all reads
                unmappedKmers += localUnmappedKmers;
                //std::cerr << "Thread " << k << " processed " << locallyProcessedReads << " reads\n"; 
            }));

            }

          }   


          // Wait for all of the threads to finish
          for ( auto& thread : threads ){ thread.join(); }

          auto end = std::chrono::steady_clock::now();
          auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
          auto nsec = sec.count();
          auto rate = (nsec > 0) ? readNum / sec.count() : 0;
          std::cerr << "\nOverall rate: " << rate << " reads / s\n";
          std::cerr << "\n" << std::endl;
          rhash.dumpCountsToFile(countsFile);

          // Total kmers
          size_t totalCount = 0;
          for (auto i : boost::irange(size_t(0), rhash.kmers().size())) {
              totalCount += rhash.atIndex(i);
          }

          bfs::path countInfoFilename(countsFile);
          countInfoFilename.replace_extension(".count_info");

          std::ofstream countInfoFile(countInfoFilename.string());
          countInfoFile << "total_reads\t" << totalCount << "\n";
          countInfoFile << "mapped\t" << totalCount - unmappedKmers << "\n";
          countInfoFile << "unmapped\t" << unmappedKmers << "\n";
          countInfoFile << "mapped_ratio\t" << 
                           (totalCount / static_cast<double>(totalCount + unmappedKmers)) << "\n";
          countInfoFile.close();

          std::cerr << "There were " << totalCount << ", kmers; " << unmappedKmers << " could not be mapped\n";
          std::cerr << "Mapped " << 
                       (totalCount / static_cast<double>(totalCount + unmappedKmers)) * 100.0 << "% of the kmers\n";
          end = std::chrono::steady_clock::now();
          sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
          nsec = sec.count();

          std::cerr << "Total counting time [including file i/o]: " << nsec << " seconds.\n";

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
