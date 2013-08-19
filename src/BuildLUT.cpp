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

#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for_each.h"
#include "tbb/task_scheduler_init.h"

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>
#include <jellyfish/parse_dna.hpp>

#include "LookUpTableUtils.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "CountDBNew.hpp"
#include "ezETAProgressBar.hpp"

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
  jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);

  std::vector<std::thread> threads;
  std::vector<TranscriptList> transcriptsForKmer;

  size_t numTranscripts = tgmap.numTranscripts();
  //std::vector<TranscriptInfo*> transcripts;
  //transcripts.resize(numTranscripts, nullptr);

  auto merLen = transcriptHash.kmerLength();
  bool done {false};
  std::atomic<size_t> numRes {0};
  std::atomic<size_t> nworking{numThreads-1};

  // Start the thread that will print the progress bar
  std::cerr << "number of kmers : " << transcriptHash.size() << "\n";
  std::cerr << "Building transcript <-> kmer lookup tables \n";
  threads.push_back( std::thread( [&numRes, &nworking, numTranscripts] () {
    size_t lastCount = numRes;
    ez::ezETAProgressBar show_progress(numTranscripts);
    show_progress.start();
    while ( nworking > 0 ) {
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

  std::mutex tmut;

  transcriptsForKmer.resize( transcriptHash.size() );
  tbb::concurrent_queue<ContainingTranscript> q;

  /**
   * spawn off a thread to build the kmer look up table
   */
  threads.push_back(std::thread(
    [&q, &transcriptsForKmer, &nworking]() {
      ContainingTranscript ct;
      while( nworking > 0 ) {
        while( q.try_pop(ct) ) {
          transcriptsForKmer[ct.kmerID].push_back(ct.transcriptID);
        }
      }
    })
  );

  transcriptsForKmer.resize( transcriptHash.size() );

  using LUTTools::TranscriptInfo;
  tbb::concurrent_queue<TranscriptInfo*> tq;
  //boost::lockfree::queue<TranscriptInfo*> tq(numTranscripts);

  /**
   * spawn off a thread to dump the transcript lookup table to file
   */
  threads.push_back(std::thread(
    [&tq, &transcriptsForKmer, &nworking, tlutfname]() {
      std::ofstream tlutstream(tlutfname, std::ios::binary);
      size_t numRec = 0;
      tlutstream.write(reinterpret_cast<const char*>(&numRec), sizeof(numRec));

      TranscriptInfo* ti = nullptr;
      while( nworking > 0 ) {
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


  // Start the desired number of threads to parse the transcripts
  // and build our data structure.
  for (size_t i = 0; i < numThreads - 1; ++i) {

    threads.push_back( std::thread( 
      [&numRes, &q, &tq, &tgmap, &parser, &transcriptHash, &nworking, 
       &transcriptIndex, &transcriptsForKmer, &tmut, merLen]() -> void {

        // Each thread gets it's own stream
        jellyfish::parse_read::thread stream = parser.new_thread();
        jellyfish::parse_read::read_t* read;
        auto INVALID = transcriptHash.INVALID;
        bool useCanonical{transcriptIndex.canonical()};

        // while there are transcripts left to process
        while ( (read = stream.next_read()) ) { 
          // The transcript name
          std::string fullHeader(read->header, read->hlen);
          std::string header = fullHeader.substr(0, fullHeader.find(' '));
          
          // The transcript sequence
          // strip the newlines from the sequence and put it into a string (newSeq)
          std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
          auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
          auto readLen = std::distance( seq.begin(), newEnd );
          std::string newSeq(seq.begin(), seq.begin() + readLen);
    
          // Lookup the ID of this transcript in our transcript -> gene map
          auto transcriptIndex = tgmap.findTranscriptID(header); 
          bool valid = ((transcriptIndex != tgmap.INVALID) and 
                        (readLen > merLen));
          auto geneIndex = tgmap.gene(transcriptIndex);

          if ( not valid ) { continue; }
          ++numRes;

          size_t numKmers {static_cast<size_t>(readLen) - merLen + 1};

          TranscriptInfo* tinfo = new TranscriptInfo;
          tinfo->name = header;
          tinfo->transcriptID = transcriptIndex;
          tinfo->geneID = geneIndex;
          tinfo->length = readLen;
          //tinfo->kmers.resize(numKmers);

          // Iterate over the kmers
          ReadLength effectiveLength(0);
          for ( auto offset : boost::irange( size_t(0), numKmers) ) { 
            // the kmer and it's uint64_t representation
            auto mer = newSeq.substr( offset, merLen );             
            auto binMer = jellyfish::parse_dna::mer_string_to_binary(mer.c_str(), merLen );
            
            if (useCanonical) {
              auto rcMer = jellyfish::parse_dna::reverse_complement(binMer, merLen);
              binMer = std::min(binMer, rcMer);
            }

            auto binMerId = transcriptHash.id(binMer);
            // Only count and track kmers which should be considered
            if ( binMerId != INVALID ) {
              auto tcount = transcriptHash.atIndex(binMerId);
              if ( tcount > 0 ) {
               ContainingTranscript c{binMerId, static_cast<TranscriptID>(transcriptIndex)};
               q.push(c);
               // for boost lockfree queue
               // while(!q.push( c ));
             }
           }
         }

         tq.push(tinfo);
         // for boost lockfree queue
         // while(!tq.push(tinfo));
         //transcripts[transcriptIndex] = tinfo;
       }
              
       --nworking;
     }) );

  }

  // Wait for all of the threads to finish
  for ( auto& thread : threads ){ thread.join(); }

  std::cerr << "writing kmer lookup table . . . ";
  std::cerr << "table size = " << transcriptsForKmer.size() << " . . . ";
  LUTTools::dumpKmerLUT(transcriptsForKmer, klutfname);
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
      ("genes,g", po::value< std::vector<string> >(), "gene sequences")
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

    std::vector<string> genesFile = vm["genes"].as<std::vector<string>>();
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
      boost::archive::binary_oarchive oa(ofs);
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

}
