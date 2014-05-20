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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "BiasIndex.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "CountDBNew.hpp"
#include "SailfishConfig.hpp"
#include "VersionChecker.hpp"
#include "LookUpTableUtils.hpp"

/**
* Type aliases
*/
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;

int main(int argc, char* argv[]) {
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  try{

   bool poisson = false;

   uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ("cfg,f", po::value< string >(), "config file")
    ;

    po::options_description config("Configuration");
    config.add_options()
      ("counts,c", po::value<string>(), "count file")
      ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
      ("bias,b", po::value<string>(), "bias index prefix (without .bin/.dict)")
      ("output,o", po::value<string>(), "output file")
      //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
      ("lutfile,l", po::value<string>(), "Lookup table prefix")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);

    if (vm.count("help")) {
	  std::cout << programOptions << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    string hashFile = vm["counts"].as<string>();
    //std::vector<string> genesFile = vm["genes"].as<std::vector<string>>();
    //string transcriptHashFile = vm["thash"].as<string>();
    string sfIndexBase = vm["index"].as<string>();
    string sfIndexFile = sfIndexBase+".sfi";
    string sfTrascriptCountFile = sfIndexBase+".sfc";
    bfs::path outputFilePath = bfs::path(vm["output"].as<string>());
    string lutprefix = vm["lutfile"].as<string>();
    auto tlutfname = lutprefix + ".tlut";
    auto klutfname = lutprefix + ".klut";

    TranscriptGeneMap tgm;
    { // read the serialized transcript <-> gene map from file
      string tgmFile = sfIndexBase+".tgm";
      std::cerr << "Reading the transcript <-> gene map from [" <<
                   tgmFile << "]\n";
      std::ifstream ifs(tgmFile, std::ios::binary);
      boost::archive::binary_iarchive ia(ifs);
      ia >> tgm;
      std::cerr << "done\n";
    }

    std::cerr << "Reading transcript index from [" << sfIndexFile << "] . . .";
    auto sfIndex = PerfectHashIndex::fromFile( sfIndexFile );
    auto del = []( PerfectHashIndex* h ) -> void { /*do nothing*/; };
    auto sfIndexPtr = std::shared_ptr<PerfectHashIndex>( &sfIndex, del );
    std::cerr << "done\n";

    /*
    std::cerr << "Reading transcript counts from [" << sfTrascriptCountFile << "] . . .";
    auto transcriptHash = CountDBNew::fromFile(sfTrascriptCountFile, sfIndexPtr);
    std::cerr << "done\n";
    */

    // the READ hash
    std::cerr << "Reading read counts from [" << hashFile << "] . . .";
    auto hash = CountDBNew::fromFile( hashFile, sfIndexPtr );
    std::cerr << "done\n";
    //const std::vector<string>& geneFiles{genesFile};
    auto merLen = sfIndex.kmerLength();

    KmerIDMap transcriptsForKmer;

    // Get the kmer look-up-table from file
    LUTTools::readKmerLUT(klutfname, transcriptsForKmer);

    // For each kmer
    size_t unique = 0;
	size_t uniqueAndMappable = 0;
	size_t mappable = 0;
	size_t totalCountedKmers = 0;
	size_t uniquelyMappedCount = 0;
    for(auto kmerIdx : boost::irange(size_t{0}, hash.size())) {
    	size_t countInReads = hash.atIndex(kmerIdx);
    	// Does is map to only a single transcript?
    	if (countInReads >= 0) {
    		unique += (transcriptsForKmer[kmerIdx].size() <= 1) ? 1 : 0;
    	}
	  	if (countInReads >= 1 ) {
    		uniqueAndMappable += (transcriptsForKmer[kmerIdx].size() <= 1) ? 1 : 0;
    		uniquelyMappedCount += (transcriptsForKmer[kmerIdx].size() <= 1) ? countInReads : 0;
    		mappable++;
    	}
    	totalCountedKmers += countInReads;
    }

    std::cerr << "There were " << unique << " unique k-mers\n";
    std::cerr << "This is " << 100.0 * (static_cast<double>(unique)/hash.size()) << "% of all kmers\n";
	std::cerr << "There were " << uniqueAndMappable << " unique & mapped k-mers\n";
    std::cerr << "This is " << 100.0 * (static_cast<double>(uniqueAndMappable)/hash.size()) << "% of all kmers\n";
	std::cerr << "This is " << 100.0 * (static_cast<double>(uniqueAndMappable)/mappable) << "% of mappable kmers\n";
	std::cerr << "By count " << 100.0 * (static_cast<double>(uniquelyMappedCount)/totalCountedKmers) << "% of kmers are uniquely mapped\n";
	size_t readLength = static_cast<int>((static_cast<double>(hash.totalLength()) / hash.numLengths()));
	std::cerr << "read length is " << readLength << "\n";
	size_t totalPossibleKmers = hash.numLengths() * (readLength - hash.kmerLength());

// mapped  8132179875
// unmapped        449494400
// mapped_ratio    0.950229
// unique  0.370477
// unique_and_mapped       0.178816

    std::ofstream outputFile(outputFilePath.string(), std::ofstream::out | std::ofstream::app);
	outputFile << "mapped\t" << totalCountedKmers << "\n";
	outputFile << "unmapped\t" << totalPossibleKmers - totalCountedKmers << "\n";
	outputFile << "mapped_ratio\t" << static_cast<double>(totalCountedKmers) / totalPossibleKmers << "\n";
    outputFile << "unique\t" << unique << "\n";
	outputFile << "unique_and_mapped\t" << uniqueAndMappable << "\n";
	outputFile << "uniquely_mapped_count\t" << uniquelyMappedCount << "\n";
	outputFile << "mappable\t" << mappable << "\n";
	outputFile << "number_of_kmers\t" << hash.size() << "\n";
	outputFile << "kmer_len\t" << hash.kmerLength() << "\n";
    outputFile.close();

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " " << "analyze_kmer_tradeoff" << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " " << "analyze_kmer_tradeoff" << " --help\nExiting.\n";
    std::exit(1);
  }

    return 0;
}
