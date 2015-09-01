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
#include <future>
#include <list>
#include <unordered_map>

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

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

#include "ReadLibrary.hpp"

#include "CountDBNew.hpp"
#include "cmph.h"

#include "PerfectHashIndex.hpp"
#include "LibraryFormat.hpp"

enum class MerDirection : std::int8_t { FORWARD = 1, REVERSE = 2, BOTH = 3 };

template <typename ParserT>
bool countKmers(std::vector<std::string> readFiles,
                //ParserT& parser,
                PerfectHashIndex& phi, CountDBNew& rhash, size_t merLen,
                bool discardPolyA, ReadStrandedness direction, std::atomic<uint64_t>& numReadsProcessed,
                std::atomic<uint64_t>&unmappedKmers, std::atomic<uint64_t>& readNum, size_t numThreads) {

  using std::string;
  using std::cerr;
  using std::vector;
  using std::thread;
  using std::atomic;
  using stream_manager = jellyfish::stream_manager<char**>;
  using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;


  char** fnames = new char*[readFiles.size()];// fnames[1];
  for (size_t i = 0; i < readFiles.size(); ++i) {
      fnames[i] = const_cast<char*>(readFiles[i].c_str());
  }
  //fnames[0] = const_cast<char*>(readFile.c_str());

  // Create a jellyfish parser
  const int concurrentFile = std::min(readFiles.size(), numThreads);
  stream_manager streams(fnames, fnames + readFiles.size(), concurrentFile);

  size_t maxReadGroupSize{100};
  sequence_parser parser(4*numThreads, maxReadGroupSize, concurrentFile, streams);

  boost::timer::auto_cpu_timer t(cerr);
  auto start = std::chrono::steady_clock::now();
  bool canonical = phi.canonical();

  atomic<size_t> fileReadNum{0};
  vector<thread> threads;
  atomic<bool> notDone{true};
  // Start the desired number of threads to parse the reads
  // and build our data structure.
  for (size_t k = 0; k < numThreads; ++k) {
    size_t threadIdx = k;


    threads.emplace_back(thread(
            [&parser, &readNum, &fileReadNum, &rhash, &start, &phi, &unmappedKmers, &k, discardPolyA, threadIdx, direction, merLen]() mutable -> bool {
                    using BinMer = uint64_t;
                    vector<BinMer> fwdMers;
                    vector<BinMer> revMers;

                    BinMer lshift{2 * (merLen - 1)};
                    BinMer masq{(1UL << (2 * merLen)) - 1};
                    BinMer cmlen;

                    size_t numKmers = 0;
                    size_t numRemaining = 0;
                    size_t fCount = 0; size_t rCount = 0;
                    auto dir = direction;

                    auto INVALID = phi.INVALID;

                    jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t> polyA(merLen);
                    polyA.polyA();

                    jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t> polyT(merLen);
                    polyT.polyT();

                    jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t>  kmer(merLen);
                    jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t>  rkmer(merLen);

                    uint64_t localUnmappedKmers{0};
                    uint64_t locallyProcessedReads{0};

                    // while there are transcripts left to process
                    while (true) {
                        sequence_parser::job j(parser);
                        // If this job is empty, then we're done
                        if (j.is_empty()) {
                            unmappedKmers += localUnmappedKmers;
                            return true ;
                        }

                        for (size_t i=0; i < j->nb_filled; ++i) {
                            ++readNum; ++locallyProcessedReads; ++fileReadNum;
                            if (readNum % 250000 == 0) {
                                auto end = std::chrono::steady_clock::now();
                                auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
                                auto nsec = sec.count();
                                auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                                cerr << "processed " << readNum << " reads (" << rate << ") reads/s\r\r";
                            }

                            const char* start     = j->data[i].seq.c_str();
                            uint32_t readLen      = j->data[i].seq.size();
                            const char* const end = start + readLen;

                            // reset all of the counts
                            fCount = rCount = numKmers = 0;
                            cmlen = 0;
                            kmer.polyA(); rkmer.polyA();
                            dir = direction;

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
                                char base = *start; ++start;
                                auto c = jellyfish::mer_dna::code(base);
                                kmer.shift_left(c);
                                rkmer.shift_right(jellyfish::mer_dna::complement(c));
                                switch(c) {
                                    case jellyfish::mer_dna::CODE_IGNORE: break;
                                    case jellyfish::mer_dna::CODE_COMMENT:
                                           std::cerr << "ERROR: unexpected character " << base << " in read!\n";

                                    // Fall through
                                    case jellyfish::mer_dna::CODE_RESET:
                                           cmlen = 0;
                                           kmer.polyA(); rkmer.polyA();
                                           break;

                                default:
                                  // count if the kmer is valid in the forward and
                                  // reverse directions
                                  if(++cmlen >= merLen) {
                                    cmlen = merLen;

                                    if (discardPolyA and (kmer == polyA or rkmer == polyA)) {
                                      ++numKmers; --numRemaining;
                                      break;
                                    }

                                    // dispatch on the direction
                                    switch (dir) {
                                       // We're certain that more kmers map in the forward direction
                                       // so we only consider the rest of the read in this direction.
                                       case ReadStrandedness::S:
                                        // get the index of the forward kmer
                                        binMerId = phi.index(kmer.get_bits(0, 2*merLen));
                                        if (binMerId != INVALID) {
                                          rhash.incAtIndex(binMerId);
                                          ++fCount;
                                        }
                                        ++numKmers; --numRemaining;
                                        break;
                                       // end case FORWARD

                                       // We're certain that more kmers map in the reverse direction
                                       // so we only consider the rest of the read in this direction.
                                       case ReadStrandedness::A:
                                          // get the index of the forward kmer
                                          rMerId = phi.index(rkmer.get_bits(0, 2*merLen));
                                          if (rMerId != INVALID) {
                                            rhash.incAtIndex(rMerId);
                                            ++rCount;
                                          }
                                          ++numKmers; --numRemaining;
                                          break;
                                       // end case REVERSE

                                       case ReadStrandedness::U:
                                          // form the new kmer and it's reverse complement

                                          // Find the index of the forward kmer and determine
                                          // whether or not to count it.
                                          binMerId = phi.index(kmer.get_bits(0, 2*merLen));
                                          fwdMers[fCount] = binMerId;
                                          fCount += (binMerId != INVALID);

                                          // Find the index of the reverse kmer and determine
                                          // whether or not to count it.
                                          rMerId = phi.index(rkmer.get_bits(0, 2*merLen));
                                          revMers[rCount] = rMerId;
                                          rCount += (rMerId != INVALID);

                                          ++numKmers; --numRemaining;

                                          // Determine if we need to continue looking at both directions
                                          dir = (fCount > (rCount + numRemaining)) ? ReadStrandedness::S :
                                                (rCount > (fCount + numRemaining)) ? ReadStrandedness::A : ReadStrandedness::U;

                                          switch (dir) {
                                            case ReadStrandedness::S:
                                              for (auto i : boost::irange(size_t(0), fCount)) { rhash.incAtIndex(fwdMers[i]);
                                              }
                                              break;
                                            case ReadStrandedness::A:
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
                          case ReadStrandedness::U:
                            if (dir == ReadStrandedness::U) {
                              for (auto i : boost::irange(size_t(0), fCount)) { rhash.incAtIndex(fwdMers[i]);
                              }
                            }
                            count = fCount;
                            break;

                          // More things mapped in the forward direction
                          case ReadStrandedness::S:
                            count = fCount; break;

                          // More things mapped in the reverse direction
                          case ReadStrandedness::A:
                            count = rCount; break;
                        }

                        // the number of unmapped kmers is just the total kmers in this read
                        // minus the number that mapped.
                        localUnmappedKmers += (numKmers - count);

                        //producer.finishedWithRead(s);
                    } // end job
                } // end parse all reads
            }));

         }

          // Wait for all of the threads to finish
          for ( auto& thread : threads ){ thread.join(); }
          cerr << "\n";
          return true;
}

/**
* source: http://stackoverflow.com/questions/10890242/get-the-status-of-a-stdfuture
*/
template<typename R>
bool futureIsReady(std::future<R>& f) {
      return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

//int mainCount( int argc, char *argv[] ) {
int mainCount( uint32_t numThreads,
               const std::string& sfIndexBase,
               const std::vector<ReadLibrary>& readLibraries,
               const std::string& countsFile,
               bool discardPolyA) {


    using std::vector;
    using std::map;
    using std::string;
    using std::cerr;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    /*
    uint32_t maxThreads = std::thread::hardware_concurrency();
    std::vector<string> undirReadFiles;// = vm["reads"].as<std::vector<string>>();
    std::vector<string> fwdReadFiles;// = vm["forward"].as<std::vector<string>>();
    std::vector<string> revReadFiles;// = vm["reverse"].as<std::vector<string>>();

    po::options_description generic("Sailfish count options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
    ("reads,r", po::value<std::vector<string>>(&undirReadFiles)->multitoken(), "List of files containing \"undirected\" reads")
    ("reverse,R", po::value<std::vector<string>>(&revReadFiles)->multitoken(), "List of files containing \"anti-sense\" reads")
    ("forward,F", po::value<std::vector<string>>(&fwdReadFiles)->multitoken(), "List of files containing \"sense\" reads")
    ("counts,c", po::value<string>(), "File where Sailfish read count is written")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
    ("polya,a", po::bool_switch(), "polyA/polyT k-mers should be discarded")
    ;

    po::variables_map vm;
    */

    try {
        size_t numActors = numThreads;
        string sfTrascriptIndexFile = sfIndexBase+".sfi";

        /*
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
        bool discardPolyA = vm["polya"].as<bool>();
        size_t numActors = vm["threads"].as<uint32_t>();

        // auto cpuset = new cpu_set_t;
        // CPU_ZERO(cpuset);
        // CPU_SET(0, cpuset);
        // if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), cpuset)) {
        //   std::cerr << "COULD NOT SET PROCESSOR AFFINITY!!!\n";
        //   std::exit(1);
        // }
        */

        std::cerr << "reading index . . . ";
        auto phi = PerfectHashIndex::fromFile(sfTrascriptIndexFile);
        std::cerr << "done\n";
        std::cerr << "index contained " << phi.numKeys() << " kmers\n";

        size_t nkeys = phi.numKeys();
        size_t merLen = phi.kmerLength();

        tbb::task_scheduler_init init(numActors);
        std::vector<std::thread> threads;

        auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
        auto phiPtr = std::shared_ptr<PerfectHashIndex>(&phi, del);

        CountDBNew rhash( phiPtr );

        std::atomic<uint64_t> readNum{0};
        std::atomic<uint64_t> processedReads{0};
        std::atomic<uint64_t> numReadsProcessed{0};
        std::atomic<uint64_t> unmappedKmers{0};

        { // create a scope --- the timer will be destructed at the end

          boost::timer::auto_cpu_timer t(std::cerr);
          auto start = std::chrono::steady_clock::now();
          std::vector<std::tuple<const std::string&, ReadStrandedness, CountDBNew*>> filesToProcess;

          struct StrandHasher {
              size_t operator()(const ReadStrandedness& s) const {
                  switch(s) {
                      case ReadStrandedness::SA:
                          return 0;
                      case ReadStrandedness::AS:
                          return 1;
                      case ReadStrandedness::A:
                          return 2;
                      case ReadStrandedness::S:
                          return 3;
                      case ReadStrandedness::U:
                          return 4;
                      default:
                          return 4;
                  }

              }
          };

          std::unordered_map<ReadStrandedness, std::vector<std::string>, StrandHasher> strandFileMap;

          for (auto& rl : readLibraries) {
              auto& libFmt = rl.format();
              auto& mate1ReadFiles = rl.mates1();
              auto& mate2ReadFiles = rl.mates2();
              auto& unmatedReadFiles = rl.unmated();
              ReadStrandedness orientation;
              if (libFmt.type == ReadType::PAIRED_END) {
                  for (auto& readFile : mate1ReadFiles) {
                      switch (libFmt.strandedness) {
                      case ReadStrandedness::SA:
                      case ReadStrandedness::S:
                          orientation = ReadStrandedness::S;
                          break;
                      case ReadStrandedness::AS:
                      case ReadStrandedness::A:
                          orientation = ReadStrandedness::A;
                          break;
                      case ReadStrandedness::U:
                          orientation = ReadStrandedness::U;
                          break;
                      }
                      filesToProcess.push_back(make_tuple(std::ref(readFile), orientation, &rhash));
                      strandFileMap[orientation].push_back(readFile);
                  }

                  for (auto& readFile : mate2ReadFiles) {
                      switch (libFmt.strandedness) {
                      case ReadStrandedness::AS:
                      case ReadStrandedness::S:
                          orientation = ReadStrandedness::S;
                          break;
                      case ReadStrandedness::SA:
                      case ReadStrandedness::A:
                          orientation = ReadStrandedness::A;
                          break;
                      case ReadStrandedness::U:
                          orientation = ReadStrandedness::U;
                          break;
                      }
                      filesToProcess.push_back(make_tuple(std::ref(readFile), orientation, &rhash));
                      strandFileMap[orientation].push_back(readFile);
                  }
              } else if (libFmt.type == ReadType::SINGLE_END) {

                  for (auto& readFile : unmatedReadFiles) {
                      auto orientation = libFmt.strandedness;
                      filesToProcess.push_back(make_tuple(std::ref(readFile), orientation, &rhash));
                      strandFileMap[orientation].push_back(readFile);
                  }
              }
          }

        using stream_manager = jellyfish::stream_manager<char**>;
        using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;

        size_t freeThreads = numActors;
        size_t maxThreadsPerFile = 8;
        size_t numProcessed = 0;
        bool done{false};
        std::list<std::tuple<std::future<bool>, size_t>> workingFiles;
/*
        do {
            if (freeThreads > 0 and numProcessed < filesToProcess.size() ) {
                size_t cthreads = std::min(maxThreadsPerFile, freeThreads);

                freeThreads -= cthreads;

                auto& countJob = filesToProcess[numProcessed];
                auto& readFile = std::get<0>(countJob);
                auto orientation = std::get<1>(countJob);
                auto& countHash = *std::get<2>(countJob);
                cerr << "file " << readFile << ": ";
                cerr << "cthreads = " << cthreads << "\n";
                auto myFuture = std::async(std::launch::async,
                        countKmers<sequence_parser>, readFile, //std::ref(parser),
                        std::ref(phi), std::ref(countHash), merLen, discardPolyA,
                        orientation, std::ref(numReadsProcessed),
                        std::ref(unmappedKmers), std::ref(readNum), cthreads);

                workingFiles.push_back(make_tuple(std::move(myFuture), cthreads));
                ++numProcessed;
            } else {
                std::chrono::seconds oneSec(1);
                std::this_thread::sleep_for(oneSec);
                // Go through the list of files currently being processed, and remove any that
                // are finished
                workingFiles.remove_if([&freeThreads](std::tuple<std::future<bool>, size_t>& e) {
                        bool isDone = futureIsReady(std::get<0>(e));
                        if (isDone) { freeThreads += std::get<1>(e); }
                        return isDone;
                        });
            }
        } while(!workingFiles.empty() or numProcessed < filesToProcess.size());
*/
        size_t totFiles{filesToProcess.size()};
        std::vector<std::future<bool>> futures;
        for (auto& kv : strandFileMap) {
            size_t curFiles = kv.second.size();
            if (curFiles > 0) {
                size_t cthreads = std::ceil( numActors * (static_cast<float>(curFiles) / totFiles));
                std::cerr << "processing files [";
                for (auto& fn : kv.second) {
                    std::cerr << " " << fn;
                }
                std::cerr << " ] using " << cthreads << " threads\n";

                auto myFuture = std::async(std::launch::async,
                        [&] () -> bool {
                        return countKmers<sequence_parser>(
                        kv.second,
                        std::ref(phi), std::ref(rhash), merLen, discardPolyA,
                        kv.first, std::ref(numReadsProcessed),
                        std::ref(unmappedKmers), std::ref(readNum), cthreads); });
                futures.emplace_back(std::move(myFuture));
            }
        }
        for (auto& f : futures) { f.get(); }

        /*
          for (auto& countJob : filesToProcess) {
              auto& readFile = std::get<0>(countJob);
              auto orientation = std::get<1>(countJob);
              auto& countHash = *std::get<2>(countJob);
              cerr << "file " << readFile << ": \n";

              namespace bfs = boost::filesystem;
              bfs::path filePath(readFile);

              char** fnames = new char*[1];// fnames[1];
              fnames[0] = const_cast<char*>(readFile.c_str());

              // Create a jellyfish parser
              const int concurrentFile{1};

              using stream_manager = jellyfish::stream_manager<char**>;
              using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;

              stream_manager streams(fnames, fnames + 1, concurrentFile);

              size_t maxReadGroupSize{100};
              sequence_parser parser(4*numActors, maxReadGroupSize, concurrentFile, streams);
*/
        /*

        std::vector<std::string> files;
        for (auto& countJob : filesToProcess) {
            files.push_back(std::get<0>(countJob));
        }
        auto orientation = std::get<1>(filesToProcess[0]);
        auto& countHash = *std::get<2>(filesToProcess[0]);
              countKmers<sequence_parser>(//parser,
                      files,
                      phi, countHash, merLen, discardPolyA,
                      orientation , numReadsProcessed,
                      unmappedKmers, readNum, numActors);
         */
/*
              cerr << "\n";
          }
            */
          auto end = std::chrono::steady_clock::now();
          auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
          auto nsec = sec.count();
          auto rate = (nsec > 0) ? readNum / sec.count() : 0;
          std::cerr << "\nOverall rate: " << rate << " reads / s\n";
          std::cerr << "\n" << std::endl;
          rhash.dumpCountsToFile(countsFile);

          // Total kmers
          size_t mappedKmers= 0;
          for (auto i : boost::irange(size_t(0), rhash.kmers().size())) {
            mappedKmers += rhash.atIndex(i);
          }

          bfs::path countInfoFilename(countsFile);
          countInfoFilename.replace_extension(".count_info");

          size_t totalCount = mappedKmers + unmappedKmers;
          std::ofstream countInfoFile(countInfoFilename.string());
          countInfoFile << "total_kmers\t" << totalCount << "\n";
          countInfoFile << "mapped\t" << mappedKmers << "\n";
          countInfoFile << "unmapped\t" << unmappedKmers << "\n";
          countInfoFile << "mapped_ratio\t" <<
                           (mappedKmers / static_cast<double>(totalCount)) << "\n";
          countInfoFile.close();

          std::cerr << "There were " << totalCount << ", kmers; " << unmappedKmers << " could not be mapped\n";
          std::cerr << "Mapped " <<
                       (mappedKmers / static_cast<double>(totalCount)) * 100.0 << "% of the kmers\n";
          end = std::chrono::steady_clock::now();
          sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
          nsec = sec.count();

          std::cerr << "Total counting time [including file i/o]: " << nsec << " seconds.\n";

        }

    } catch (po::error &e) {
        std::cerr << "Program Options Error : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception &e) {
        std::cerr << "ERROR: sailfish count (subordinate command) invoked improperly.\n";
        std::exit(1);
    }

    return 0;
}
