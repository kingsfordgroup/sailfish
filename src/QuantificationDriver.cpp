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
#include <cstdlib>
#include <fstream>
#include <vector>
#include <thread>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

using std::string;

int mainCount(int argc, char* argv[]);
int runIterativeOptimizer(int argc, char* argv[]);

int runKmerCounter(const std::string& sfCommand, 
                   uint32_t numThreads,
                   const std::string& indexBase, 
                   const std::vector<string>& readFiles, 
                   const std::string& countFileOut) {

    std::stringstream argStream;
    argStream << sfCommand << " ";
    argStream << "--index " << indexBase << " ";
    argStream << "--counts " << countFileOut << " ";
    argStream << "--threads " << numThreads << " ";

    argStream << "--reads ";
    for (auto& rfile : readFiles) {
        argStream << rfile << " ";
    }

    std::string argString = argStream.str();
    boost::trim(argString);

    // Run the Sailfish counting phase as an separate process.
    // This will force the mmapped memory to be cleaned up.
    auto pid = fork();
    if (pid == 0) { //child
        // Create the argv array for the main call to Sailfish
        std::vector<std::string> argStrings;
        boost::split(argStrings, argString, boost::is_any_of(" "));

        char ** args = new char*[argStrings.size()];
        for (size_t i : boost::irange({0}, argStrings.size())) {
            // args[i] = new char[argStrings[i].length()+1];
            // std::copy(argStrings[i].begin(), argStrings[i].end(), args[i]);
            // args[i][argStrings[i].length()] = '\0';
            args[i] = const_cast<char*>(argStrings[i].c_str());
        }

        std::cerr << "In Sailfish Counting thread. Counting read kmers\n";
        int ret = mainCount(argStrings.size(), args);
        delete [] args;
        std::exit(ret);

    } else if (pid < 0) { // fork failed!
        
        std::cerr << "FATAL ERROR: Failed to spawn Sailfish process. Exiting\n";
        std::exit(-1);
    } else { // parent

        int status = -1;
        waitpid(pid, &status, 0); // wait on the Jellyfish process
        std::cerr << "Sailfish terminated with return code " << status << "\n";
        assert( status == 0 );
    }

    return 0;
}

    
int runSailfishEstimation(const std::string& sfCommand, 
                          uint32_t numThreads, 
                          const boost::filesystem::path& countFile, 
                          const boost::filesystem::path& indexBasePath,
                          //const std::string& tgmap, 
                          size_t iterations, 
                          const boost::filesystem::path& lookupTableBase, 
                          const boost::filesystem::path& outFilePath,
                          bool noBiasCorrect) {

    std::stringstream argStream;
    argStream << sfCommand << " ";
    // argStream << "estimate ";
    if (noBiasCorrect) {
        argStream << "--no_bias_correct ";
    }
    argStream << "--index " << indexBasePath.string() << " ";
    argStream << "--counts " << countFile.string() << " ";
    argStream << "--threads " << numThreads << " ";
    //argStream << "--tgmap " << tgmap << " ";
    argStream << "--lutfile " << lookupTableBase.string() << " ";
    argStream << "--iterations " << iterations << " ";
    argStream << "--out " << outFilePath.string();

    std::string argString = argStream.str();
    boost::trim(argString);

    int ret = 0;
    auto estThread = std::thread(
        [&]() -> void {
            std::vector<std::string> argStrings;
            boost::split(argStrings, argString, boost::is_any_of(" "));

            char ** args = new char*[argStrings.size()];
            for (size_t i : boost::irange({0}, argStrings.size())) {
                args[i] = const_cast<char*>(argStrings[i].c_str());
            }

            std::cerr << "In Sailfish estimation thread. Estimating transcript abundances\n";
            int ret = runIterativeOptimizer(argStrings.size(), args);
            delete [] args;
        }
    );

    /*
    std::cerr << "cmdline: " << argString << "\n";
    auto sfEstFile = popen(argString.c_str(), "r");
    int estRet = pclose(sfEstFile);
    std::cerr << "Sailfish estimate finished with status: " << estRet << "\n";
    return estRet;
    */
   estThread.join();
   return ret;
}

int mainQuantify( int argc, char *argv[] ) {

    using std::string;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    string sfCommand = argv[0];
    uint32_t maxThreads = std::thread::hardware_concurrency();
    bool noBiasCorrect = false;

    po::options_description generic("Sailfish quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "Sailfish index [output of the \"Sailfish index\" command")        
    ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")    
    //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("out,o", po::value<string>(), "Basename of file where estimates are written")
    ("iterations,n", po::value<size_t>()->default_value(30), "number of iterations to run the optimzation")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
    ("force,f", po::bool_switch(), "Force the counting phase to rerun, even if a count databse exists." )
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

        bfs::path indexBasePath(vm["index"].as<string>());
        bfs::path outputBasePath(vm["out"].as<string>());
        //string tgmap = vm["tgmap"].as<string>();
        //string tgmap = indexBase+".tgm";
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        std::vector<string> readFiles = vm["reads"].as<std::vector<string>>();
        bool force = vm["force"].as<bool>();

        /*
        ("index,i", po::value<string>(), "transcript index file [Sailfish format]")
        ("reads,r", po::value<std::vector<string>>()->multitoken(), "List of files containing reads")
        ("counts,c", po::value<string>(), "File where Sailfish read count is written")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
        */

        if (bfs::exists(outputBasePath) and !bfs::is_directory(outputBasePath)) {
            std::cerr << "The provided output path [" << outputBasePath << "] " <<
                         "already exists and is not a directory\n.";
            std::cerr << "Please either provide a different path or " <<
                         "delete the existing file.\n";
            std::exit(1);
        }

        bool mustRecount = false;
        if (!bfs::exists(outputBasePath)) {
            mustRecount = true;
            try {
                bool success = bfs::create_directory(outputBasePath);
                if (!success) { throw std::runtime_error("unspecified error creating file."); }
            } catch ( std::exception& e ) {
                std::cerr << "Creation of " << outputBasePath << " failed [" << e.what() << "]\n.";
                std::cerr << "Exiting.\n";
                std::exit(1);
            }
        }

        //string countFileOut = outputBase + ".sfc";
        bfs::path countFilePath(outputBasePath); countFilePath /= "reads.sfc";
        bfs::path indexPath(indexBasePath); indexPath /= "transcriptome";

        mustRecount = (force or !boost::filesystem::exists(countFilePath));
        if (mustRecount) {
            runKmerCounter(sfCommand, numThreads, indexPath.string(), readFiles, countFilePath.string());
        }

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

        bfs::path lutBasePath(indexBasePath); lutBasePath /= "transcriptome";
        bfs::path estFilePath(outputBasePath); estFilePath /= "quant.sf";
        size_t iterations = vm["iterations"].as<size_t>();
        runSailfishEstimation(sfCommand, numThreads, countFilePath, indexPath,
                              iterations, lutBasePath, estFilePath, noBiasCorrect);

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
