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

#include "ReadProducer.hpp"

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
#include "StreamingSequenceParser.hpp"

enum class MerDirection : std::int8_t { FORWARD = 1, REVERSE = 2, BOTH = 3 };

template <typename ParserT>
bool countKmers(ParserT& parser, PerfectHashIndex& phi, CountDBNew& rhash, size_t merLen,
                bool discardPolyA, MerDirection direction, std::atomic<uint64_t>& numReadsProcessed,
                std::atomic<uint64_t>&unmappedKmers, std::atomic<uint64_t>& readNum, size_t numThreads) {

  using std::string;
  using std::cerr;
  using std::vector;
  using std::thread;
  using std::atomic;

  boost::timer::auto_cpu_timer t(cerr);
  auto start = std::chrono::steady_clock::now();
  bool canonical = phi.canonical();

  atomic<size_t> fileReadNum{0};
  vector<thread> threads;
  // Start the desired number of threads to parse the reads
  // and build our data structure.
  for (size_t k = 0; k < numThreads; ++k) {
    size_t threadIdx = k;


    threads.emplace_back(thread(
            [&parser, &readNum, &fileReadNum, &rhash, &start, &phi, &unmappedKmers, &k, discardPolyA, threadIdx, direction, merLen]() mutable -> void {
                    using BinMer = uint64_t;
                    vector<BinMer> fwdMers;
                    vector<BinMer> revMers;

                    BinMer lshift{2 * (merLen - 1)};
                    BinMer masq{(1UL << (2 * merLen)) - 1};
                    BinMer cmlen, kmer, rkmer;

                    size_t numKmers = 0;
                    size_t numRemaining = 0;
                    size_t fCount = 0; size_t rCount = 0;
                    auto dir = direction;

                    auto INVALID = phi.INVALID;
                    char* as = new char[merLen];
                    for (auto i : boost::irange(size_t{0}, merLen)) { as[i] = 'A'; }

                    auto polyA = jellyfish::parse_dna::mer_string_to_binary(as, merLen );

                    uint64_t localUnmappedKmers{0};
                    uint64_t locallyProcessedReads{0};

                    ReadProducer<ParserT> producer(parser);

                    ReadSeq* s;

                    while (producer.nextRead(s)) {
                        ++readNum; ++locallyProcessedReads; ++fileReadNum;
                        if (readNum % 250000 == 0) {
                            auto end = std::chrono::steady_clock::now();
                            auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start);
                            auto nsec = sec.count();
                            auto rate = (nsec > 0) ? fileReadNum / sec.count() : 0;
                            cerr << "processed " << readNum << " reads (" << rate << ") reads/s\r\r";
                        }
                        const char* start     = s->seq;
                        uint32_t readLen      = s->len;
                        const char* const end = s->seq + readLen;

                        // reset all of the counts
                        fCount = rCount = numKmers = 0;
                        cmlen = kmer = rkmer = 0;
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

                                    if (discardPolyA and (kmer == polyA or rkmer == polyA)) {
                                      ++numKmers; --numRemaining;
                                      break;
                                    }

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

                        producer.finishedWithRead(s);

                } // end parse all reads
                unmappedKmers += localUnmappedKmers;
                delete [] as;
            }));

         }

          // Wait for all of the threads to finish
          for ( auto& thread : threads ){ thread.join(); }
          cerr << "\n";
}


int mainCount( int argc, char *argv[] ) {

    using std::vector;
    using std::map;
    using std::string;
    using std::cerr;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

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
    ("reverse,R", po::value<std::vector<string>>(&fwdReadFiles)->multitoken(), "List of files containing \"sense\" reads")
    ("forward,F", po::value<std::vector<string>>(&revReadFiles)->multitoken(), "List of files containing \"anti-sense\" reads")
    ("counts,c", po::value<string>(), "File where Sailfish read count is written")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
    ("polya,a", po::bool_switch(), "polyA/polyT k-mers should be discarded")
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
        bool discardPolyA = vm["polya"].as<bool>();

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

        CountDBNew rhash( phiPtr );

        std::atomic<uint64_t> readNum{0};
        std::atomic<uint64_t> processedReads{0};
        std::atomic<uint64_t> numReadsProcessed{0};
        std::atomic<uint64_t> unmappedKmers{0};

        map<MerDirection, vector<string>> fnameMap{
            {MerDirection::BOTH, undirReadFiles},
            {MerDirection::FORWARD, fwdReadFiles},
            {MerDirection::REVERSE, revReadFiles}};

        { // create a scope --- the timer will be destructed at the end

          boost::timer::auto_cpu_timer t(std::cerr);
          auto start = std::chrono::steady_clock::now();

        for (auto& df : fnameMap) {

          auto& direction = df.first;

          cerr << "Parsing ";
          switch (direction) {
              case MerDirection::BOTH:
                cerr << "undirected ";
                break;
              case MerDirection::FORWARD:
                cerr << "sense ";
                break;
              case MerDirection::REVERSE:
                cerr << "antisense ";
                break;
          }
          cerr << "reads\n";

          for (auto& readFile : df.second) {
            cerr << "file " << readFile << ": \n";

            namespace bfs = boost::filesystem;
            bfs::path filePath(readFile);

            // If this is a regular file, then use the Jellyfish parser
            if (bfs::is_regular_file(filePath)) {

              char** fnames = new char*[1];// fnames[1];
              fnames[0] = const_cast<char*>(readFile.c_str());

              jellyfish::parse_read parser(fnames, fnames+1, 5000);

              countKmers<jellyfish::parse_read>(
                         parser, phi, rhash, merLen, discardPolyA,
                         direction, numReadsProcessed,
                         unmappedKmers, readNum, numActors);

            } else { // If this is a named pipe, then use the kseq-based parser
              vector<bfs::path> paths{readFile};
              StreamingReadParser parser(paths);
              parser.start();

              countKmers<StreamingReadParser>(
                         parser, phi, rhash, merLen, discardPolyA,
                         direction, numReadsProcessed,
                         unmappedKmers, readNum, numActors);
            }

            cerr << "\n";
          }

        }

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
          countInfoFile << "total_reads\t" << totalCount << "\n";
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
        std::cerr << "ERROR: " << argv[0] << " count invoked improperly.\n";
        std::cerr << "Usage\n";
        std::cerr << "=====\n";
        std::cout << generic << std::endl;
        std::exit(1);
    }

    return 0;
}
