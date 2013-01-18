
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

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/find.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include "g2logworker.h"
#include "g2log.h"

#include "boost/lockfree/fifo.hpp"
#include "threadpool.hpp"

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "tclap/CmdLine.h"
#include "affymetrix_utils.hpp"
#include "kfitf.hpp"
#include "linear_system_builder.hpp"
#include "transcript_segmenter.hpp"
#include "set_multicover_solver.hpp"
#include "iterative_optimizer.hpp"
#include "utils.hpp"
#include "genomic_feature.hpp"
#include "CountDB.hpp"

using affymetrix::utils::AffyEntry;

using seqan::DnaString;


typedef std::tuple< std::string, double, double > record;


struct IdProbT{
  size_t id;
  double prob;
};

struct SeqProbT{
  std::string seq;
  double prob;
};

std::ostream& operator<< ( std::ostream& os, const std::vector<IdProbT>& p ) {
  for ( auto& e : p ) {
    os << "\t" << e.id << ", " << e.prob;
  }
  return os;
}

typedef std::unordered_map<std::string, double> NaiveKmerMapT;

NaiveKmerMapT readSeqFile ( const std::string& fname, uint32_t& merSize ) {
  NaiveKmerMapT ents;
  std::ifstream affyFile(fname);

  auto update = [&]( const std::string& k, double v ) -> void {
   auto kIt = ents.find(k);
   // If the k-mer already exists in the map
   if ( kIt != ents.end() ) { kIt->second += v; } else { ents[k] = v; }
  };

  std::string seq;
  double expression, bg;
  std::string line;
  size_t i = 0;
  size_t failCount = 0;

  while ( std::getline(affyFile, line) ) {
    std::stringstream ls(line);
    ls >> seq >> expression >> bg ;
    if ( i % 10000 == 0 ) { std::cerr << "i = " << i << "\n"; }
    if ( !ls.fail() ) {
      auto bIndex = 0;
      auto eIndex = seq.length() - (merSize + 1);

      while (bIndex <= eIndex) {
        update(seq.substr(bIndex, merSize), (expression / bg) );
        ++bIndex;
      }
    } else {
      failCount += 1;
    }

    ++i;
  }

  std::cerr << failCount << " lines could not be correctly parsed and were ignored\n";

  affyFile.close();
  return ents;
}

bool dumpMap ( const std::string& fname, const NaiveKmerMapT& kmerMap ) {

  std::ofstream ofile(fname);

  for ( const auto& kv : kmerMap ) {
    ofile << kv.first << "\t" << kv.second << "\n";
  }

  ofile.close();
}

int command1(int argc, char* argv[] ) {
  using std::string;
  try{
    TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    TCLAP::ValueArg<string> affyFileName("a", "affy", "affymetrix probes with intensities", true, "", "string");
    TCLAP::ValueArg<string> outputFileName("o", "output", "file to which output should be written", true, "", "string");
    TCLAP::ValueArg<uint32_t> merSize("m", "mersize", "size of k-mers to count", true, 10, "uint32_t");
    cmd.add(affyFileName);
    cmd.add(outputFileName);
    cmd.add(merSize);

    cmd.parse(argc, argv);
    auto entries = readSeqFile( affyFileName.getValue(), merSize.getValue() );
    std::cerr << "found " << entries.size() << " distinct kmers \n";
    dumpMap ( outputFileName.getValue(), entries );
  } catch (TCLAP::ArgException &e){
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
}

template <typename KeyT, typename ValueT>
std::unordered_map<KeyT, ValueT> readIntoHash( const std::string& fname, bool reverse = false ) {
  std::unordered_map<KeyT, ValueT> m;
  std::ifstream ifile(fname.c_str());

  KeyT k{};
  ValueT v{};

  std::unordered_map< char, char > rc{ {'A','T'}, {'T','A'}, {'C','G'}, {'G','C'} };

  auto reverseComplement = [&]( const std::string& s ) {
    std::string rs( s.size(), ' ' );
    size_t i = 0;
    for ( const auto& c : s ) {
      rs[i] = rc[c];
      ++i;
    }
    return rs;
  };

  while ( ifile >> k >> v ) {
    k = reverse ? reverseComplement(k) : k;
    m[k] = v;
  }

  ifile.close();
  return m;
}


/** Uses C style file I/O because it really does seem to be
    that much faster! **/
std::vector<SeqProbT> readSeqs( const std::string& fname ) {
  std::vector<SeqProbT> s;

  std::ifstream ifile(fname.c_str(), std::ios_base::binary);

  uint64_t numLines=0;
  uint32_t probeLength=0;

  ifile.read( reinterpret_cast<char*>(&numLines), sizeof(numLines));
  ifile.read( reinterpret_cast<char*>(&probeLength), sizeof(probeLength));

  char cseq[255];// = new char[255];//std::string seq;
  float expression, bg;

  s.reserve(numLines);


  //size_t offset = + sizeof(numLines) + sizeof(probeLength);
  //offset += (numLines-100)*((sizeof(char)*probeLength) + sizeof(expression) + sizeof(bg));
  //ifile.seekg( offset );
  //numLines = 100;

  for (size_t i = 0; i < numLines; ++i ) {
    uint32_t si = 0;
    //ifile.read( reinterpret_cast<char*>(&si), sizeof(si) );
    //std::cout << "size = " << si << "\n";
    ifile.read( cseq, sizeof(char)*probeLength);
    std::string seq(cseq, probeLength);
    ifile.read( reinterpret_cast<char*>(&expression), sizeof(expression) );
    ifile.read( reinterpret_cast<char*>(&bg), sizeof(bg) );
    s.push_back( SeqProbT{seq, (expression/bg)} );
    assert( seq[probeLength-1] != ' ' );
    std::cout << "> " << i << ":" << expression << ":" << bg << "\n";
    std::cout << seq << std::endl;

    if ( i % 100000 == 0 ) { std::cerr << i << "\n"; }
  }

  ifile.close();
  return s;
}


int command2(int argc, char* argv[] ) {
  using std::string;
  try{
    TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    TCLAP::ValueArg<string> seqFileName("s", "seq", "affymetrix sequences", true, "", "string");
    TCLAP::ValueArg<string> affyFileName("a", "affy", "affymetrix kmers", true, "", "string");
    TCLAP::ValueArg<string> rnaFileName("r", "rna", "rna-seq kmers", true, "", "string");

    cmd.add(seqFileName);
    cmd.add(affyFileName);
    cmd.add(rnaFileName);


    cmd.parse(argc, argv);
    std::cerr << "Reading seqs\n";
    auto seqs = readSeqs(seqFileName.getValue());
    std::cerr << "Reading affy kmers\n";
    auto amap = readIntoHash<std::string, double>(affyFileName.getValue());
    std::cerr << "Reading rna kmers\n";
    auto rmap = readIntoHash<std::string, uint32_t>(rnaFileName.getValue(), true);

    std::unordered_map<std::string, std::tuple<uint32_t,double>> omap;
    double X(0), Y(0), EX(0), EY(0), EXY(0), EX2(0), EY2(0);
    size_t N(0);

    for ( auto& kv : rmap ){
      if ( amap.find(kv.first) != amap.end() ) {
        ++N;
        X = kv.second; Y = amap[kv.first];
        EX += X; EY += Y;

        for ( auto& s : seqs ) {
          if ( s.seq.find(kv.first) != std::string::npos ) {
            double percent = s.prob / kv.second;
          }
        }

        omap[kv.first] = std::make_tuple( kv.second, amap[kv.first] );
        std::cerr << "omap[" << kv.first << "] = (" << kv.second << "," << amap[kv.first]  << ")\n";
      }
    }

    EX /= N; EY /= N;
    for ( const auto& kv : omap ) {
      X = std::get<0>(kv.second) - EX;
      Y = std::get<1>(kv.second) - EY;
      EXY += X*Y; EX2 += X*X; EY2 += Y*Y;
    }
    std::cerr << "N = " << N << ", EXY = " << EXY << ", EX = " << EX << ", EY = " << EY << ", EX2 = " << EX2 << ", EY2 = " << EY2 << "\n";
    double pearsonR = (EXY) / (sqrt(EX2*EY2)+1e-20);
    std::cerr << "pearson's r = " << pearsonR << "\n";
  } catch (TCLAP::ArgException &e){
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
}

using std::string;
using std::shared_ptr;
//using namespace cppa;

struct KmerMappedInfoT {
    string name;
    float mappedCount;
    float transcriptLength;
};

template <typename hash_t, typename stream_t>
bool countKmersInTranscripts( 
  hash_t& hash, 
  hash_t& transcriptHash, 
  std::unordered_set<std::string>& uniqueMers, 
  KFITFCalculator& KFITFCalc, 
  stream_t& inputStream, 
  const std::string& ofname 
  ) {

  size_t merLen = hash.get_mer_len();

  size_t ctr{0};
  std::unordered_map<std::string, KmerMappedInfoT> geneIDs;

  
  std::cerr << "mer length is " << merLen << "\n";
  size_t queryCounter = 0;

  size_t numActors = 26;

  struct TranscriptEntry {
    std::string header;
    std::string seq;
  };

  boost::threadpool::pool tp(numActors);
  boost::lockfree::fifo< KmerMappedInfoT* > q;
  boost::lockfree::fifo< TranscriptEntry* > workQueue;

  std::mutex resLock;
  size_t k = 0;

  auto update = [&geneIDs] ( KmerMappedInfoT* km ) -> bool {
    auto pos = geneIDs.find(km->name);
    if ( pos == geneIDs.end() ) {
      geneIDs[km->name] = KmerMappedInfoT{km->name, km->mappedCount, km->transcriptLength};
    } else {

      pos->second.mappedCount += km->mappedCount;
      pos->second.transcriptLength += km->transcriptLength;
      /*
      if ( km->mappedCount > pos->second.mappedCount ) {
        pos->second = KmerMappedInfoT{km->name, km->mappedCount, km->transcriptLength};
      }
      */
    }
  };

  auto readResults = [&q, &tp, &update, &resLock] () -> bool {
    size_t numRes = 0;
    while (tp.active() + tp.pending() > 1) {
      KmerMappedInfoT* km;
      while( q.dequeue(km) ) {

        update(km);

        delete km;
        ++numRes;
        std::cerr << "#res = " << numRes << "\n";
      }
      boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
    }
    return true;
  };

  tp.schedule(readResults);

  bool done = false;

  for (size_t i = 0; i < numActors-1; ++i ) {
  auto task = [&hash, &q, &workQueue, &uniqueMers, &transcriptHash, &done, merLen]() -> void {
    TranscriptEntry* te;
    while ( !done ) {
      while (workQueue.dequeue(te) ) {
      string seq = te->seq;
      auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
      auto readLen = std::distance( seq.begin(), newEnd );

      string header = te->header;
      //auto splitPos = header.find('.');
      auto geneName = header;//.substr(0, splitPos);
      size_t offset = 0;
      size_t mapped = 0;
      float weight = 0;
      std::vector< std::tuple<uint64_t, uint32_t> > mappedUniquemers;
      std::unordered_map<string, size_t> tfmap;
      size_t maxTermFreq{0};

      // For each kmer in the read
      while ( offset < readLen - merLen + 1 )  {
        // Check if this kmer is unique 
        auto mer = seq.substr( offset, merLen );
        auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
        /*
        if ( std::binary_search( uniqueMers.begin(), uniqueMers.end(), binMer ) ) { 
          // If so, then it will help determine this transcript's abundance
          // mapped += hash[ mer.c_str()];
          mappedUniquemers.push_back( std::make_tuple(binMer, hash[mer.c_str()]) );
        }
        */
        //bool hasUniquemer = uniqueMers.find(mer) != uniqueMers.end();
        if ( transcriptHash[mer.c_str()] < 15 ) {
          tfmap[mer] += 1;
          maxTermFreq = (maxTermFreq > tfmap[mer] ) ? maxTermFreq : tfmap[mer];
          //if ( transcriptHash[mer.c_str()] == 1 ) { uniqueMapped += hash[mer]; }
        }
        ++offset;
      }

      /*
      float inv = 1.0 / mappedUniquemers.size();
      for ( auto& uniquemer : mappedUniquemers ) {
        weight += inv * std::get<1>(uniquemer);
      }
      float weightLen = merLen * mappedUniquemers.size();
      */      

      
      float invMaxTermFreq = 1.0 / maxTermFreq;
      offset = 0;
      float weightLen{0.0f};
      for ( auto& mc : tfmap ) {
        auto invWeight = (1.0f / (1.0 + transcriptHash[ mc.first.c_str() ]));
        weight += mc.second * hash[ mc.first.c_str() ] *  invWeight;//KFITFCalc.itf(mc.first);
        weightLen += mc.second * invWeight;// KFITFCalc.itf(mer);
      }
      /*
      while ( offset < readLen - merLen )  {
        auto mer = seq.substr( offset, merLen );
        weight += hash[mer.c_str()] * KFITFCalc.itf(mer);
        weightLen += KFITFCalc.itf(mer);
        ++offset;
      }
      */

      KmerMappedInfoT* km = new KmerMappedInfoT{geneName, weight, weightLen};
      q.enqueue(km);
      }
      boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
      }
    };

    tp.schedule(task);
  }

  size_t numReads = 0;
  jellyfish::read_parser::read_t *read;
  while ( (read = inputStream.next_read()) ) {//}&& queryCounter < 10000 ) {
    std::string oseq( read->seq_s, std::distance(read->seq_s, read->seq_e)-1 );
    //auto newEnd  = std::remove( readStr.begin(), readStr.end(), '\n' );
    //auto readLen = std::distance( readStr.begin(), newEnd );
    std::string header( read->header, read->hlen );
    auto te = new TranscriptEntry{ header, oseq };
    workQueue.enqueue( te );
    ++numReads;
   }
   done = true;

    // auto splitPos = header.find('|');
    // auto geneName = header.substr(1, splitPos-1);



  tp.wait();

  std::cerr << "Database contained " << ctr << " entries comprising " << geneIDs.size() << " distinct genes\n";

  std::ofstream ofile(ofname);

  for ( const auto& kv : geneIDs ) {
    ofile << kv.first << '\t' << kv.second.transcriptLength << '\t' << kv.second.mappedCount << '\n';
  }

  ofile.close();


  //cppa::await_all_others_done();
  return true;
}

std::unordered_set<string> loadUniqueSet( const std::string& fname ) {

  std::ifstream ifile(fname.c_str());
  size_t d = 0;
  size_t numMers = 0;
  ifile.read( reinterpret_cast<char*>(&d), sizeof(d) );
  ifile.read( reinterpret_cast<char*>(&numMers), sizeof(numMers) );
  std::cerr << "Using " << d << " mers";
  char str[d];

  std::unordered_set<std::string> us;
  us.reserve( numMers );
  size_t i = 0;
  while ( !ifile.eof() ) {
    ifile.read( str, sizeof(str) );
    us.insert(str);
    ++i;
    if ( i % 10000 == 0 ) { std::cerr << i << "\n"; }
  }
  ifile.close();
  return us;

}


int weightedCount(int argc, char* argv[] ) {
  using std::string;
  try{
    TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    TCLAP::ValueArg<string> genesFile("g", "genes", "gene sequences", true, "", "string");
    TCLAP::ValueArg<string> transcriptHashFile("t", "thash", "transcript jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> hashFile("j", "jfhash", "jellyfish hash file", true, "", "string");
    //TCLAP::ValueArg<string> uniqueFile("u", "unique", "file containing unique kmers", true, "", "string");
    TCLAP::ValueArg<string> outFile("o", "output", "file to store output in", true, "", "string");
    //TCLAP::ValueArg<string> affyFileName("a", "affy", "affymetrix probes with intensities", true, "", "string");
    //TCLAP::ValueArg<string> outputFileName("o", "output", "file to which output should be written", true, "", "string");
    //TCLAP::ValueArg<uint32_t> merSize("m", "mersize", "size of k-mers to count", true, 10, "uint32_t");
    cmd.add(genesFile);
    cmd.add(hashFile);
    // cmd.add(uniqueFile);
    cmd.add(transcriptHashFile);
    cmd.add(outFile);
    //cmd.add(affyFileName);
    //cmd.add(outputFileName);
    //cmd.add(merSize);

    cmd.parse(argc, argv);

    std::cerr << "reading file " << genesFile.getValue() << "\n";
    const char* fnames[] = { genesFile.getValue().c_str() };

    KFITFCalculator kf(18);
    /*
    {
      jellyfish::parse_read parser( fnames, fnames+1, 100);
      jellyfish::parse_read::thread stream = parser.new_thread();

      std::cerr << "reading file " << hashFile.getValue() << "\n";

      kf.computeITF(stream);
    }
    */

    jellyfish::parse_read parser( fnames, fnames+1, 100);
    jellyfish::parse_read::thread stream = parser.new_thread();

    auto uniqueKmers = std::unordered_set<string>();//loadUniqueSet( uniqueFile.getValue() );
    //auto uniqueKmers = loadUniqueSet( uniqueFile.getValue() );
    mapped_file transcriptDB( transcriptHashFile.getValue().c_str() );
    transcriptDB.random().will_need();
    char typeTrans[8];
    memcpy(typeTrans, transcriptDB.base(), sizeof(typeTrans));

    mapped_file dbf( hashFile.getValue().c_str() );
    dbf.random().will_need();
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));

    for ( size_t i = 0; i < sizeof(type); ++i ) {
      assert(type[i] == typeTrans[i]);
    }

    if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
      raw_inv_hash_query_t hash(dbf);
      raw_inv_hash_query_t transcriptHash(transcriptDB);
      countKmersInTranscripts( hash, transcriptHash, uniqueKmers, kf, stream, outFile.getValue() );
    } else if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type))) {

      hash_query_t hash(hashFile.getValue().c_str());
      hash_query_t transcriptHash(transcriptHashFile.getValue().c_str());
      countKmersInTranscripts( hash, transcriptHash, uniqueKmers, kf, stream, outFile.getValue() );
   }

  } catch (TCLAP::ArgException &e){
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
}


int runNNLS(int argc, char* argv[] ) {
  using std::string;
  try{
    TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    TCLAP::ValueArg<string> genesFile("g", "genes", "gene sequences", true, "", "string");
    TCLAP::ValueArg<string> hashFile("j", "jfhash", "jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> transcriptHashFile("t", "thash", "transcript jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> transcriptGeneMap("m", "tgmap", "file that maps transcripts to genes", true, "", "string");
    //TCLAP::ValueArg<string> uniqueFile("u", "unique", "file containing unique kmers", true, "", "string");
    //TCLAP::ValueArg<string> matFileName("m", "matout", "file to store A in", true, "", "string");
    //TCLAP::ValueArg<string> bFileName("b", "bout", "file to store b in", true, "", "string");
    cmd.add(transcriptHashFile);
    cmd.add(genesFile);
    cmd.add(hashFile);
    cmd.add(transcriptGeneMap);
    //cmd.add(uniqueFile);
    //cmd.add(matFileName);
    //cmd.add(bFileName);

    cmd.parse(argc, argv);

    typedef GenomicFeature<TranscriptGeneID> CustomGenomicFeature;

    std::ifstream transcriptGeneFile(transcriptGeneMap.getValue() );
    std::vector< CustomGenomicFeature > features;
    CustomGenomicFeature feat;
    std::cerr << "parsing gtf file [" << transcriptGeneMap.getValue()  << "] . . . ";
    while ( transcriptGeneFile >> feat ) { 
      features.push_back(feat);
    };
    std::cerr << "done\n";
    transcriptGeneFile.close();
    std::cerr << "building transcript to gene map . . .";
    auto tgm = utils::transcriptToGeneMapFromFeatures( features );
    std::cerr << "done\n";
    /*std::ifstream transcriptGeneFile( transcriptGeneMap.getValue() );
    auto tgm = utils::readTranscriptToGeneMap( transcriptGeneFile );
    transcriptGeneFile.close();
    */
    
    mapped_file transcriptDB( transcriptHashFile.getValue().c_str() );
    transcriptDB.random().will_need();
    char typeTrans[8];
    memcpy(typeTrans, transcriptDB.base(), sizeof(typeTrans));
    
    mapped_file dbf( hashFile.getValue().c_str() );
    dbf.random().will_need();
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));

    if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
      std::cerr << "Raw hash not supported!\n";
      std::abort();
      /*
      raw_inv_hash_query_t hash(dbf);
      raw_inv_hash_query_t transcriptHash(transcriptDB);
      TranscriptSegmenter<raw_inv_hash_query_t> ts( genesFile.getValue(), hash, tgm );
      ts.findOverlappingSegments( hash.get_mer_len(), transcriptHash );
      
      LinearSystemBuilder<raw_inv_hash_query_t> lsb( genesFile.getValue(), hash);
      lsb.writeToFiles( matFileName.getValue(), bFileName.getValue() );
      */
    } else if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type))) {
      hash_query_t hash(hashFile.getValue().c_str());
      hash_query_t transcriptHash(transcriptHashFile.getValue().c_str());
      TranscriptSegmenter<hash_query_t> ts( genesFile.getValue(), hash, tgm );
      ts.findOverlappingSegments( hash.get_mer_len(), transcriptHash );
      /*
      LinearSystemBuilder<hash_query_t> lsb( genesFile.getValue(), hash);
      lsb.writeToFiles( matFileName.getValue(), bFileName.getValue() );
      */
   }

  } catch (TCLAP::ArgException &e){
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
}

int runSetCover(int argc, char* argv[] ) {
  using std::string;
  try{
    TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    TCLAP::ValueArg<string> genesFile("g", "genes", "gene sequences", true, "", "string");
    TCLAP::ValueArg<string> hashFile("j", "jfhash", "jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> transcriptHashFile("t", "thash", "transcript jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> outputFile("o", "output", "output file", true, "", "string");
    cmd.add(transcriptHashFile);
    cmd.add(genesFile);
    cmd.add(hashFile);
    cmd.add(outputFile);

    cmd.parse(argc, argv);

    mapped_file transcriptDB( transcriptHashFile.getValue().c_str() );
    transcriptDB.random().will_need();
    char typeTrans[8];
    memcpy(typeTrans, transcriptDB.base(), sizeof(typeTrans));

    mapped_file dbf( hashFile.getValue().c_str() );
    dbf.random().will_need();
    char type[8];
    memcpy(type, dbf.base(), sizeof(type));

   if(!strncmp(type, jellyfish::raw_hash::file_type, sizeof(type))) {
      raw_inv_hash_query_t hash(dbf);
      raw_inv_hash_query_t transcriptHash(transcriptDB);
      //std::cerr << "transcriptHash size is " << transcriptHash.get_distinct() << "\n";
      SetMulticoverSolver<raw_inv_hash_query_t, raw_inv_hash_query_t> solver( hash, transcriptHash );
      solver( genesFile.getValue(), outputFile.getValue() );
    } else if(!strncmp(type, jellyfish::compacted_hash::file_type, sizeof(type))) {
      hash_query_t hash(hashFile.getValue().c_str());
      hash_query_t transcriptHash(transcriptHashFile.getValue().c_str());
      std::cerr << "transcriptHash size is " << transcriptHash.get_distinct() << "\n";
      SetMulticoverSolver<hash_query_t, hash_query_t> solver( hash, transcriptHash );
      solver( genesFile.getValue(), outputFile.getValue() );
   }

  } catch (TCLAP::ArgException &e){
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }
}

int runIterativeOptimizer(int argc, char* argv[] ) {
  using std::string;
  namespace po = boost::program_options;

  try{

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ("cfg,f", po::value< string >(), "config file")
    ;

    po::options_description config("Configuration");
    config.add_options()
      ("genes,g", po::value< std::vector<string> >(), "gene sequences")
      ("count,c", po::value<string>(), "count file")
      ("thash,t", po::value<string>(), "transcript jellyfish hash file")
      ("output,o", po::value<string>(), "output file")
      ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
      ("iterations,i", po::value<size_t>(), "number of iterations to run the optimzation")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);

    if ( vm.count("help") ){
      std::cout << "Sailfish\n";
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

    /*
    TCLAP::CmdLine cmd("Sailfish", ' ', "1.0");
    TCLAP::MultiArg<string> genesFile("g", "genes", "gene sequences", true, "string");
    //TCLAP::ValueArg<string> hashFile("j", "jfhash", "jellyfish hash file", true, "", "string");

    TCLAP::ValueArg<string> orderFile("s", "orderFile", "file storing order of the kmers", true, "", "string");
    TCLAP::ValueArg<string> hashFile("c", "count", "count file", true, "", "string");

    TCLAP::ValueArg<string> transcriptHashFile("t", "thash", "transcript jellyfish hash file", true, "", "string");
    TCLAP::ValueArg<string> outputFile("o", "output", "output file", true, "", "string");
    TCLAP::ValueArg<string> transcriptGeneMap("m", "tgmap", "file that maps transcripts to genes", true, "", "string");
    TCLAP::ValueArg<size_t> numIter("i", "iterations", "number of iterations to run the optimization", false, 3, "size_t");
    cmd.add(transcriptHashFile);
    cmd.add(genesFile);
    cmd.add(transcriptGeneMap);
    cmd.add(orderFile);
    cmd.add(hashFile);
    cmd.add(outputFile);
    cmd.add(numIter);

    cmd.parse(argc, argv);
    */

    string transcriptGeneMap = vm["tgmap"].as<string>();
    string hashFile = vm["count"].as<string>();
    std::vector<string> genesFile = vm["genes"].as<std::vector<string>>();
    string transcriptHashFile = vm["thash"].as<string>();
    string outputFile = vm["output"].as<string>();
    size_t numIter = vm["iterations"].as<size_t>();

    typedef GenomicFeature<TranscriptGeneID> CustomGenomicFeature;

    std::ifstream transcriptGeneFile(transcriptGeneMap );
    std::vector< CustomGenomicFeature > features;
    CustomGenomicFeature feat;
    std::cerr << "parsing gtf file [" << transcriptGeneMap  << "] . . . ";
    while ( transcriptGeneFile >> feat ) { 
      features.push_back(feat);
    };
    std::cerr << "done\n";
    transcriptGeneFile.close();
    std::cerr << "building transcript to gene map . . .";
    auto tgm = utils::transcriptToGeneMapFromFeatures( features );
    std::cerr << "done\n";

    std::cerr << "Reading transcript index from [" << transcriptHashFile << "] . . .";
    auto transcriptHash = CountDB::fromFile( transcriptHashFile );
    std::cerr << "done\n";

    // the READ hash
    std::cerr << "Reading indexed read counts from [" << hashFile << "] . . .";
    CountDB hash( hashFile, transcriptHash.indexKmers(), transcriptHash.kmerLength() );    
    std::cerr << "done\n";

    const std::vector<string>& geneFiles{genesFile};
    auto merLen = transcriptHash.kmerLength();
    
    std::cerr << "Creating optimizer . . .";
    IterativeOptimizer<CountDB, CountDB> solver( hash, transcriptHash, tgm );
    std::cerr << "done\n";

    std::cerr << "optimizing for " << numIter << " iterations";
    solver.optimize( geneFiles, outputFile, numIter );

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  }

}

int main( int argc, char* argv[] ) {

  g2LogWorker logger(argv[0], "./" );
  g2::initializeLogging(&logger);
  std::cerr << "** log file being written to " << logger.logFileName().get() << "** \n";
  
  using std::string;

  try {

    std::unordered_map<std::string, std::function<int(int, char*[])>> cmds({ 
      {"cmd1", command1}, 
      {"cmd2", command2}, 
      {"wcount", weightedCount}, 
      {"nnls", runNNLS},
      {"itopt", runIterativeOptimizer},
      {"setcover", runSetCover }});

    // TCLAP::CmdLine cmd("RNASeq!", ' ', "1.0");
    // TCLAP::SwitchArg cmd1("x", "xx", "cmd1", false);
    // TCLAP::SwitchArg cmd2("y", "yy", "cmd2", false);

    // std::vector<TCLAP::Arg*> cmds;
    // cmds.push_back(&cmd1);
    // cmds.push_back(&cmd2);

    // cmd.xorAdd( cmds );

    // cmd.parse(argc, argv);
    char** argv2 = new char*[argc-1];
    argv2[0] = argv[0];
    std::copy_n( &argv[2], argc-2, &argv2[1] );
    cmds[ argv[1] ](argc-1, argv2);
    delete[] argv2;

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
  }


return 0;
}
