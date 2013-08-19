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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "BiasIndex.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "CountDBNew.hpp"
#include "CollapsedIterativeOptimizer.hpp"
#include "SailfishConfig.hpp"
#include "VersionChecker.hpp"
//#include "iterative_optimizer.hpp"
//#include "tclap/CmdLine.h"

int performBiasCorrection(boost::filesystem::path featPath, 
                          boost::filesystem::path expPath, 
                          boost::filesystem::path outPath,
                          size_t numThreads);

int runIterativeOptimizer(int argc, char* argv[] ) {
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  string cmdString = "estimate";

  try{

   bool poisson = false;
   bool noBiasCorrect = false;

   uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ("cfg,f", po::value< string >(), "config file")
    ;

    po::options_description config("Configuration");
    config.add_options()
      //("genes,g", po::value< std::vector<string> >(), "gene sequences")
      ("counts,c", po::value<string>(), "count file")
      ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
      ("bias,b", po::value<string>(), "bias index prefix (without .bin/.dict)")
      //("thash,t", po::value<string>(), "transcript jellyfish hash file")
      ("output,o", po::value<string>(), "output file")
      ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
      //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
      ("filter,f", po::value<double>()->default_value(0.0), "during iterative optimization, remove transcripts with a mean less than filter")
      ("iterations,n", po::value<size_t>(), "number of iterations to run the optimzation")
      ("lutfile,l", po::value<string>(), "Lookup table prefix")
      ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);
 
    //bool poisson = ( vm.count("poisson") ) ? true : false;
    
    if ( vm.count("version") ) {
      std::cout << "version : " << Sailfish::version <<"\n";
      std::exit(0);
    }

    if ( vm.count("help") ){
      std::cout << "Sailfish\n";
      std::cout << programOptions << std::endl;
      std::exit(0);
    }

    if ( vm.count("cfg") ) {
      std::cerr << "have detected configuration file\n";
      string cfgFile = vm["cfg"].as<string>();
      std::cerr << "cfgFile : [" << cfgFile << "]\n";
      po::store(po::parse_config_file<char>(cfgFile.c_str(), programOptions, true), vm);
    }
    po::notify(vm);

    bool computeBiasCorrection = !noBiasCorrect;
    uint32_t numThreads = vm["threads"].as<uint32_t>();
    tbb::task_scheduler_init init(numThreads);

    string hashFile = vm["counts"].as<string>();
    //std::vector<string> genesFile = vm["genes"].as<std::vector<string>>();
    //string transcriptHashFile = vm["thash"].as<string>();
    
    string sfIndexBase = vm["index"].as<string>();
    
    bfs::path sfIndexBasePath(vm["index"].as<string>());

    string sfIndexFile = sfIndexBase+".sfi";
    string sfTrascriptCountFile = sfIndexBase+".sfc";
    bfs::path outputFilePath = bfs::path(vm["output"].as<string>());

    double minMean = vm["filter"].as<double>();
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
    
    BiasIndex bidx = vm.count("bias") ? BiasIndex( vm["bias"].as<string>() ) : BiasIndex();

    std::cerr << "Creating optimizer . . .";
    CollapsedIterativeOptimizer<CountDBNew> solver(hash, tgm, bidx, numThreads);
    // IterativeOptimizer<CountDBNew, CountDBNew> solver( hash, transcriptHash, tgm, bidx );
    std::cerr << "done\n";

    if ( poisson ) {
      std::cerr << "optimizing using Poisson model\n";
      // for IterativeOptimizer
      // solver.optimizePoisson( geneFiles, outputFile );
    } else {
      size_t numIter = vm["iterations"].as<size_t>();
      std::cerr << "optimizing using iterative optimization [" << numIter << "] iterations";
      // for CollapsedIterativeOptimizer (EM algorithm)

      solver.optimize(klutfname, tlutfname, numIter, minMean );

      std::stringstream headerLines;
      headerLines << "# [sailfish version]\t" << Sailfish::version << "\n";
      headerLines << "# [kmer length]\t" << sfIndex.kmerLength() << "\n";
      headerLines << "# [using canonical kmers]\t" << (sfIndex.canonical() ? "true" : "false") << "\n";
      headerLines << "# [command]\t";
      for (size_t i : boost::irange(size_t(0), static_cast<size_t>(argc))) { headerLines << argv[i] << " "; }
      headerLines << "\n";

      solver.writeAbundances(outputFilePath, headerLines.str());
      
     
      if (computeBiasCorrection) {
        sfIndexBasePath.remove_filename();
        outputFilePath.remove_filename();        
        auto biasFeatPath = sfIndexBasePath / "bias_feats.txt";
        auto expressionFilePath = outputFilePath / "quant.sf";
        auto biasCorrectedFile = outputFilePath / "quant_bias_corrected.sf";
        std::cerr << "biasFeatPath = " << biasFeatPath << "\n";
        std::cerr << "expressionFilePath = " << expressionFilePath << "\n";
        std::cerr << "biasCorrectedFile = " << biasCorrectedFile << "\n";                
        performBiasCorrection(biasFeatPath, expressionFilePath, biasCorrectedFile, numThreads);
      }

      // for LASSO Iterative Optimizer
      //solver.optimizeNNLASSO(klutfname, tlutfname, outputFile, numIter, minMean );
      // for IterativeOptimizer
      // solver.optimize( geneFiles, outputFile, numIter, minMean );
    }

    

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " " << cmdString << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " " << cmdString << " --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}

int help(int argc, char* argv[]) {
  auto helpmsg = R"(
  ===============

  Please invoke sailfish with one of the following commands {index, quant, sf}.  
  For more inforation on the options for theses particular methods, use the -h 
  flag along with the method name.  For example:

  Sailfish index -h

  will give you detailed help information about the index command.
  )"; 

  std::cerr << "  Sailfish v" << Sailfish::version << helpmsg << "\n";
  return 1;
}

/**
 * Bonus! 
 */
int mainSailfish(int argc, char* argv[]) {

  std::cerr << R"(
   _____       _ _______      __  
  / ___/____ _(_) / __(_)____/ /_ 
  \__ \/ __ `/ / / /_/ / ___/ __ \
 ___/ / /_/ / / / __/ (__  ) / / /
/____/\__,_/_/_/_/ /_/____/_/ /_/ 
)";

  return 0;

}

//int indexMain( int argc, char* argv[] );
int mainIndex(int argc, char* argv[]);
int mainCount(int argc, char* argv[]);
int mainQuantify(int argc, char* argv[]);
int mainBuildLUT(int argc, char* argv[] );

int main( int argc, char* argv[] ) {  
  using std::string;
  namespace po = boost::program_options;


  try {

    po::options_description hidden("hidden");
    hidden.add_options()
    ("command", po::value<string>(), "command to run {index, estimate, sf}");
    
    po::options_description sfopts("Allowed Options");
    sfopts.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ;

    po::options_description all("Allowed Options");
    all.add(sfopts).add(hidden);

    std::string versionMessage = getVersionMessage();
    std::cerr << versionMessage;

    // po::options_description sfopts("Command");
    // sfopts.add_options()

    po::positional_options_description pd;
    pd.add("command", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(2, argv).options(all).positional(pd).allow_unregistered().run(), vm);

    if (vm.count("version")) {
      std::cerr << "version : " << Sailfish::version << "\n";
      std::exit(0);
    }

    if (!vm.count("command") and vm.count("help")) {
      std::cout << sfopts << std::endl;
      help(argc, argv);
      std::exit(0);
    }

    po::notify(vm);
      
    std::unordered_map<string, std::function<int(int, char*[])>> cmds({
      {"estimate", runIterativeOptimizer},
      {"index", mainIndex},
      {"buildlut", mainBuildLUT},
      {"quant", mainQuantify},
      {"count", mainCount},
      {"sf", mainSailfish}
    });

    //string cmd = argv[1];
    string cmd = vm["command"].as<string>();

    char** argv2 = new char*[argc-1];
    argv2[0] = argv[0];
    std::copy_n( &argv[2], argc-2, &argv2[1] );

    auto cmdMain = cmds.find(cmd);
    if (cmdMain == cmds.end()) {
      help( argc-1, argv2 );
    } else {
      cmdMain->second(argc-1, argv2);
    }
    delete[] argv2;

  } catch (po::error &e) {
    std::cerr << "Program Option Error (main) : [" << e.what() << "].\n Exiting.\n";
    std::exit(1);
  } catch (...) {
    std::cerr << argv[0] << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " --help\nExiting.\n";
  }

  return 0;
}
