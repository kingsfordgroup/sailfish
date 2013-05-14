#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <thread>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

using std::string;

int runKmerCounter(const std::string& sfCommand, 
                   uint32_t numThreads,
                   const std::string& indexBase, 
                   const std::vector<string>& readFiles, 
                   const std::string& countFileOut) {

    std::stringstream argStream;
    argStream << sfCommand << " ";
    argStream << "count ";
    argStream << "--index " << indexBase << " ";
    argStream << "--counts " << countFileOut << " ";
    argStream << "--threads " << numThreads << " ";

    argStream << "--reads ";
    for (auto& rfile : readFiles) {
        argStream << rfile << " ";
    }

    std::string argString = argStream.str();
    boost::trim(argString);

    std::cerr << "cmdline: " << argString << "\n";
    auto sfCountFile = popen(argString.c_str(), "r");
    int countRet = pclose(sfCountFile);
    std::cerr << "Sailfish count finished with status: " << countRet << "\n";
    return countRet;
}

    
int runSailfishEstimation(const std::string& sfCommand, 
                          uint32_t numThreads, 
                          const std::string& countFile, 
                          const std::string& indexBase,
                          const std::string& tgmap, 
                          size_t iterations, 
                          const std::string& lookupTable, 
                          const std::string& outFileBase) {

    std::stringstream argStream;
    argStream << sfCommand << " ";
    argStream << "estimate ";
    argStream << "--index " << indexBase << " ";
    argStream << "--counts " << countFile << " ";
    argStream << "--threads " << numThreads << " ";
    argStream << "--tgmap " << tgmap << " ";
    argStream << "--lutfile " << lookupTable << " ";
    argStream << "--iterations " << iterations << " ";
    argStream << "--out " << outFileBase << ".sf";

    std::string argString = argStream.str();
    boost::trim(argString);

    std::cerr << "cmdline: " << argString << "\n";
    auto sfEstFile = popen(argString.c_str(), "r");
    int estRet = pclose(sfEstFile);
    std::cerr << "Sailfish estimate finished with status: " << estRet << "\n";
    return estRet;
}

int mainQuantify( int argc, char *argv[] ) {

    using std::string;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Sailfish quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "Sailfish index [output of the \"Sailfish index\" command")        
    ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
    ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("out,o", po::value<string>(), "Basename of file where estimates are written")
    ("iterations,i", po::value<size_t>()->default_value(30), "number of iterations to run the optimzation")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
    ;
 
    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ){
          std::cout << "Sailfish quant\n";
          std::cout << generic << std::endl;
          std::exit(0);
        }

        po::notify(vm);

        string indexBase = vm["index"].as<string>();
        string outputBase = vm["out"].as<string>();
        string tgmap = vm["tgmap"].as<string>();
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        std::vector<string> readFiles = vm["reads"].as<std::vector<string>>();

        /*
        ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
        ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
        ("counts,c", po::value<string>(), "File where Sailfish read count is written")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
        */

        string countFileOut = outputBase + ".sfc";
        string sfCommand = argv[0];
        runKmerCounter(sfCommand, numThreads, indexBase, readFiles, countFileOut);

        /*
        //("genes,g", po::value< std::vector<string> >(), "gene sequences")
        ("counts,c", po::value<string>(), "count file")
        ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
        ("bias,b", po::value<string>(), "bias index prefix (without .bin/.dict)")
        //("thash,t", po::value<string>(), "transcript jellyfish hash file")
        ("output,o", po::value<string>(), "output file")
        ("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
        ("filter,f", po::value<double>()->default_value(0.0), "during iterative optimization, remove transcripts with a mean less than filter")
        ("iterations,i", po::value<size_t>(), "number of iterations to run the optimzation")
        ("lutfile,l", po::value<string>(), "Lookup table prefix")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
        */

        size_t iterations = vm["iterations"].as<size_t>();
        runSailfishEstimation(sfCommand, numThreads, countFileOut, indexBase,
                              tgmap, iterations, indexBase, outputBase);

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }

    return 0;
}
