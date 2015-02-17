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
#include <boost/lockfree/queue.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <memory>
#include <functional>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <chrono>
#include <iomanip>

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options/parsers.hpp>

#include "cereal/archives/binary.hpp"

#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for_each.h"
#include "tbb/task_scheduler_init.h"

#include "LookUpTableUtils.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "CountDBNew.hpp"
#include "ezETAProgressBar.hpp"
#include "PartitionRefiner.hpp"
#include "spdlog/spdlog.h"

using TranscriptID = uint32_t;
using KmerID = uint64_t;
using ReadLength = int64_t;
using Length = uint32_t;
using TranscriptList = std::vector<TranscriptID>;

struct ContainingTranscript{
  KmerID kmerID;
  TranscriptID transcriptID;
};

/**
 * This function builds both a kmer => transcript and transcript => kmer
 * lookup table.
 */
int buildLUTs(
  const std::vector<std::string>& transcriptFiles, //!< File from which transcripts are read
  PerfectHashIndex& transcriptIndex,               //!< Index of transcript kmers
  CountDBNew& transcriptHash,                      //!< Count of kmers in transcripts
  TranscriptGeneMap& tgmap,                        //!< Transcript => Gene map
  const std::string& tlutfname,                    //!< Transcript lookup table filename
  const std::string& klutfname,                    //!< Kmer lookup table filename
  uint32_t numThreads                              //!< Number of threads to use in parallel
  ) {

  using LUTTools::TranscriptInfo;
  using std::vector;
  using std::atomic;
  using std::max_element;
  namespace bfs = boost::filesystem;
  using tbb::blocked_range;

  auto jointLog = spdlog::get("jointLog");
  auto fileLog = spdlog::get("fileLog");

  char** fnames = new char*[transcriptFiles.size()];
  size_t z{0};
  size_t numFnames{0};
  for ( auto& s : transcriptFiles ){
    // Ugly, yes?  But this is not as ugly as the alternatives.
    // The char*'s contained in fname are owned by the transcriptFiles
    // vector and need not be manually freed.
    fnames[numFnames] = const_cast<char*>(s.c_str());
    ++numFnames;
  }
  // Create a jellyfish parser
  const int concurrentFile{1};
  using stream_manager = jellyfish::stream_manager<char**>;
  using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;
  stream_manager streams(fnames, fnames + numFnames, concurrentFile);

  vector<bfs::path> paths{transcriptFiles[0]};

  vector<std::thread> threads;
  vector<TranscriptList> transcriptsForKmer;

  size_t numTranscripts = tgmap.numTranscripts();
  vector<TranscriptInfo*> transcripts;
  transcripts.resize(numTranscripts, nullptr);

  auto merLen = transcriptHash.kmerLength();
  bool done {false};
  atomic<size_t> numRes {0};
  atomic<size_t> nworking{(numThreads > 1) ? (numThreads - 1) : 1};

  size_t maxReadGroupSize{100};
  sequence_parser parser(4*nworking, maxReadGroupSize, concurrentFile, streams);

  // Start the thread that will print the progress bar
  std::cerr << "Number of kmers : " << transcriptHash.size() << "\n";
  std::cerr << "Parsing transcripts and building k-mer equivalence classes\n";

  threads.push_back( std::thread( [&numRes, &nworking, numTranscripts] () {
    size_t lastCount = numRes;
    ez::ezETAProgressBar show_progress(numTranscripts);
    show_progress.start();
    while ( numRes < numTranscripts ) {
      auto diff = numRes - lastCount;
      if ( diff > 0 ) {
        show_progress += static_cast<unsigned int>(diff);
        lastCount += diff;
      }
      boost::this_thread::sleep_for(boost::chrono::milliseconds(100));
    }
    show_progress += static_cast<unsigned int>(numTranscripts - lastCount);
    std::cerr << "\n";
  }) );


  PartitionRefiner refiner(transcriptHash.size());
  std::mutex refinerMutex;

  std::atomic<size_t> numInvalidKmers{0};
  // Start the desired number of threads to parse the transcripts
  // and build our data structure.
  size_t numWorkers = (numThreads > 1) ? numThreads -1 : 1;
  for (size_t i = 0; i < numWorkers; ++i) {

    threads.push_back( std::thread(
      [&numRes, &tgmap, &parser, &transcriptHash, &nworking, &transcripts,
       &fileLog, &jointLog, &transcriptIndex, &transcriptsForKmer, &refiner,
       &refinerMutex, &numInvalidKmers, merLen]() -> void {

        // Each thread gets it's own stream
        //jellyfish::parse_read::thread stream = parser.new_thread();
        //jellyfish::parse_read::read_t* read;

        auto INVALID = transcriptHash.INVALID;
        bool useCanonical{transcriptIndex.canonical()};

        // while there are transcripts left to process
        while (true) { //producer.nextRead(s)) {
          sequence_parser::job j(parser);
          // If this job is empty, then we're done
          if (j.is_empty()) { --nworking; return; }

          for (size_t i=0; i < j->nb_filled; ++i) {
            // The transcript name
            std::string fullHeader(j->data[i].header);
            std::string header = fullHeader.substr(0, fullHeader.find(' '));

            // The transcript sequence; strip the newlines from the
            // sequence and put it into a string (newSeq)
            //std::string seq(s->seq_s, std::distance(read->seq_s, read->seq_e) - 1);
            //auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
            //auto readLen = std::distance( seq.begin(), newEnd );
            auto readLen = j->data[i].seq.size();
            std::string newSeq(j->data[i].seq);

          // Lookup the ID of this transcript in our transcript -> gene map
          auto transcriptIndex = tgmap.findTranscriptID(header);
          bool valid = (transcriptIndex != tgmap.INVALID);
          auto geneIndex = tgmap.gene(transcriptIndex);

          if ( not valid ) { ++numRes; continue; }
          ++numRes;

          size_t numKmers {(readLen >= merLen) ? static_cast<size_t>(readLen) - merLen + 1 : 0};

          TranscriptInfo* tinfo = new TranscriptInfo;
          tinfo->name = header;
          tinfo->transcriptID = transcriptIndex;
          tinfo->geneID = geneIndex;
          tinfo->length = readLen;
          tinfo->kmers.resize(numKmers, INVALID);

          // Iterate over the kmers
          ReadLength effectiveLength(0);
          size_t nextKmerID{0};
          size_t locallyInvalidKmers{0};

          jellyfish::mer_dna_ns::mer_base_dynamic<uint64_t> kmer(merLen);

          size_t cmlen{0};
          size_t offset{0};
          size_t numChars{newSeq.size()};
          while (offset < numChars) {
              int c = jellyfish::mer_dna::code(newSeq[offset]);
              kmer.shift_left(c);
              if (jellyfish::mer_dna::not_dna(c)) {
                  cmlen = 0;
                  ++offset;
                  ++locallyInvalidKmers;
                  continue;
              }
              if (++cmlen >= merLen) {
                  if (useCanonical) {
                      kmer.canonicalize();
                  }

                  auto binMerId = transcriptHash.id(kmer.get_bits(0, 2*merLen));
                  // For now, we "handle" k-mers with Ns by skipping them
                  if (binMerId != INVALID) {
                      tinfo->kmers[nextKmerID++] = binMerId;
                  } else {
                      ++locallyInvalidKmers;
                  }
              }
              ++offset;
          }
          // Discount k-mers containing Ns
         tinfo->kmers.resize(nextKmerID);
         if (locallyInvalidKmers > 0) {
             numInvalidKmers += locallyInvalidKmers;
             fileLog->warn() << "Transcript [" << tinfo->name << "] contains " << locallyInvalidKmers << " unhashable k-mers";
         }

            // Partition refinement is not threadsafe.  Thus,
            // we have to use a mutex here to assure that multiple
            // threads don't try to refine the partitions at the
            // same time (there may be a more efficient approach).
            refinerMutex.lock();
            refiner.splitWith(tinfo->kmers);
            refinerMutex.unlock();

             transcripts[transcriptIndex] = tinfo;
             //producer.finishedWithRead(s);

          } // end of current job
       } // while (true)

     }) );

  }

  // Wait for all of the threads to finish
  for ( auto& thread : threads ){ thread.join(); }
  threads.clear();

  if (numInvalidKmers > 0) {
      jointLog->warn() << "total unhashable k-mers: " << numInvalidKmers
                       << "\nSome target transcripts contained \'N\'s or other characters resulting in some unhashable k-mers.\n"
                       << "This may be okay, but check the log -- in the \'logs\' directory under the index --- for details.\n\n";
  }

  // For simplicity and speed, the partition refiner is allowed to use many
  // more labels than there are partitions (each newly created partition gets
  // a new label, even when the label for one of the existing partitions could
  // be used).  Here, we "compact" the labels, assuring they are [0,#part-1].
  refiner.relabel();
  const auto& membership = refiner.partitionMembership();

  boost::filesystem::path p(tlutfname);
  p = p.parent_path();
  p /= "kmerEquivClasses.bin";

  // Dump the vector of k-mer equivalence classes to file
  LUTTools::dumpKmerEquivClasses(membership, p.string());

  /**
   *   Phase 2:
   *   Build the [k] => t map and write to file
   *   Build the t => [k] map and write to file
   **/

  // Start the thread that will print the progress bar
  std::cerr << "Building the k-mer equiv. class <=> transcript mappings\n";
  numRes = 0;
  atomic<size_t> numTranscriptsRemaining{numTranscripts};

  threads.push_back( std::thread( [&numRes, numTranscripts] () {
        size_t lastCount = numRes;
        ez::ezETAProgressBar show_progress(numTranscripts);
        show_progress.start();
        while ( numRes < numTranscripts ) {
          auto diff = numRes - lastCount;
          if ( diff > 0 ) {
            show_progress += static_cast<unsigned int>(diff);
            lastCount += diff;
          }
          boost::this_thread::sleep_for(boost::chrono::milliseconds(1010));
        }
        show_progress += static_cast<unsigned int>(numTranscripts - lastCount);
        std::cerr << "\n";
      }) );


  /**
   *  For each equivalence class, we keep a list of the transcripts in which it occurs.
   *  For each transcript, we keep a list of the equivalence classes it contains.
   **/
  size_t numEquivClasses = (*max_element(membership.cbegin(), membership.cend())) + 1;
  vector<TranscriptList> transcriptsForKmerClass(numEquivClasses);
  tbb::concurrent_queue<ContainingTranscript> q;
  tbb::concurrent_queue<TranscriptInfo*> tq;

  // Spawn off a thread to build the kmer look up table
  threads.push_back(std::thread(
                                [&q, &transcriptsForKmerClass, &numTranscriptsRemaining]() {
                                  ContainingTranscript ct;
                                  while (numTranscriptsRemaining > 0 or not q.empty()) {
                                    while (q.try_pop(ct)) {
                                      transcriptsForKmerClass[ct.kmerID].push_back(ct.transcriptID);
                                    }
                                  }
                                })
                    );




  // spawn off a thread to dump the transcript lookup table to file
  threads.push_back(std::thread(
                                [&tq, &numTranscriptsRemaining, tlutfname]() -> void {
                                  std::ofstream tlutstream(tlutfname, std::ios::binary);
                                  size_t numRec = 0;
                                  tlutstream.write(reinterpret_cast<const char*>(&numRec), sizeof(numRec));

                                  TranscriptInfo* ti = nullptr;
                                  while (numTranscriptsRemaining > 0 or not tq.empty()) {
                                    while( tq.try_pop(ti) ) {
                                      LUTTools::writeTranscriptInfo(ti, tlutstream);
                                      delete ti;
                                      ++numRec;
                                    }
                                  }

                                  // go back to the beginning of the file and write the total number
                                  // of actual records
                                  tlutstream.seekp(0);
                                  tlutstream.write(reinterpret_cast<const char*>(&numRec), sizeof(numRec));
                                  tlutstream.close();
                                })
                    );

  tbb::parallel_for(blocked_range<size_t>(0, transcripts.size()),
                    [&] (blocked_range<size_t>& trange) -> void {


                      auto INVALID = transcriptHash.INVALID;
                      // For every transcript in this thread's range
                      //std::cerr << "Transcripts.size() = " << transcripts.size() << "\n";
                      //std::cerr << "trange.begin() " << trange.begin() << ", trange.end() " << trange.end() << "\n";
                      for (auto tidx = trange.begin(); tidx != trange.end(); ++tidx) {
                        //std::cerr << "tidx = " << tidx << "\n";
                        //if (tidx >= transcripts.size()) {
                        //    std::cerr << "WARNING DANGER AHHHHHH!!!!\n";
                        //}
                        auto t = transcripts[tidx];

                        // We'll transform the list of k-mers to the list of k-mer
                        // equivalence classes
                        size_t numKmers = t->kmers.size();
                        vector<KmerID> kmerEquivClasses(numKmers, INVALID);
                        auto& kmers = t->kmers;

                        // Populate the equivalence class list, and pass each transcript
                        // off to the
                        for (size_t i = 0; i < numKmers; ++i) {
                          // k-mer => k-mer class

                          size_t kmerClass = membership[kmers[i]];
                          kmerEquivClasses[i] = kmerClass;

                          // populate k-mer class => transcript index
                          if ( kmerClass != INVALID ) {
                              ContainingTranscript c{kmerClass, static_cast<TranscriptID>(t->transcriptID)};
                              q.push(c);
                          } else {
                            std::cerr << "When building the lookup table, I thought " << kmerClass << " was an invalid k-mer class\n";
                          }

                        } // end for

                        // Now the transcript holds a vector of equivalence classes
                        std::swap(t->kmers, kmerEquivClasses);

                        tq.push(t);
                        --numTranscriptsRemaining;
                        ++numRes;
                      }

                     }
                   );

  for (auto& t : threads) { t.join(); }

  std::cerr << "writing k-mer equiv class lookup table . . . ";
  std::cerr << "table size = " << transcriptsForKmerClass.size() << " . . . ";
  LUTTools::dumpKmerLUT(transcriptsForKmerClass, klutfname);
  std::cerr << "done\n";

  return 0;
}

/**
 * This function is the main command line driver for the lookup table
 * building phase of Sailfish.  The 'buildlut' command that invokes this
 * driver is part of the 'advanced' Sailfish interface and thus, this function
 * is not usually called by the main pipeline.
 *
 * @param  argc number of tokens on command line
 * @param  argv command line tokens
 * @return      0 on success; a different value otherwise
 */
int mainBuildLUT(int argc, char* argv[] ) {
  using std::string;
  using std::vector;
  namespace po = boost::program_options;

  try{

   uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
    ;

    po::options_description config("Configuration");
    config.add_options()
      ("genes,g", po::value< vector<string> >(), "gene sequences")
      ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
      ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
      ("lutfile,l", po::value<string>(), "Lookup table prefix")
      ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);

    if ( vm.count("help") ){
      std::cout << "build_lut\n";
      std::cout << programOptions << std::endl;
      std::exit(1);
    }

    if ( vm.count("cfg") ) {
      std::cerr << "have detected configuration file\n";
      string cfgFile = vm["cfg"].as<string>();
      std::cerr << "cfgFile : [" << cfgFile << "]\n";
      po::store(po::parse_config_file<char>(cfgFile.c_str(), programOptions, true), vm);
    }
    po::notify(vm);

    uint32_t numThreads = vm["threads"].as<uint32_t>();
    tbb::task_scheduler_init init(numThreads);

    vector<string> genesFile = vm["genes"].as<vector<string>>();
    string sfIndexBase = vm["index"].as<string>();
    string sfIndexFile = sfIndexBase+".sfi";
    string sfTrascriptCountFile = sfIndexBase+".sfc";
    string lutprefix = vm["lutfile"].as<string>();
    auto tlutfname = lutprefix + ".tlut";
    auto klutfname = lutprefix + ".klut";

    TranscriptGeneMap tgmap;

    // If the user procided a GTF file, then use that to enumerate the
    // transcripts and build the transcript <-> gene map
    if (vm.count("tgmap") ) {
      string transcriptGeneMap = vm["tgmap"].as<string>();
      std::cerr << "building transcript to gene map using gtf file [" <<
                   transcriptGeneMap << "] . . .\n";
      auto features = GTFParser::readGTFFile<TranscriptGeneID>(transcriptGeneMap);
      tgmap = sailfish::utils::transcriptToGeneMapFromFeatures( features );
      std::cerr << "done\n";
    } else {
    // Otherwise, build the transcript <-> gene map directly from the
    // provided fasta file of transcripts
      std::cerr << "building transcript to gene map using transcript fasta file [" <<
                   genesFile[0] << "] . . .\n";
      tgmap = sailfish::utils::transcriptToGeneMapFromFasta(genesFile[0]);
      std::cerr << "done\n";
    }


    // save transcript <-> gene map to archive
    {
      string tgmOutFile = sfIndexBase+".tgm";
      std::cerr << "Saving transcritpt to gene map to [" << tgmOutFile << "] . . . ";
      std::ofstream ofs(tgmOutFile, std::ios::binary);
      cereal::BinaryOutputArchive oa(ofs);
      // write class instance to archive
      oa << tgmap;
      std::cerr << "done\n";
    } // archive and stream closed when destructors are called


    std::cerr << "Reading transcript index from [" << sfIndexFile << "] . . .";
    auto sfIndex = PerfectHashIndex::fromFile( sfIndexFile );
    auto del = []( PerfectHashIndex* h ) -> void { };
    auto sfIndexPtr = std::shared_ptr<PerfectHashIndex>( &sfIndex, del );
    std::cerr << "done\n";

    std::cerr << "Reading transcript counts from [" << sfTrascriptCountFile << "] . . .";
    auto transcriptHash = CountDBNew::fromFile(sfTrascriptCountFile, sfIndexPtr);
    std::cerr << "done\n";

    buildLUTs(genesFile, sfIndex, transcriptHash, tgmap, tlutfname, klutfname, numThreads);

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  }

  return 0;
}
