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

#include "LibraryFormat.hpp"

using std::string;

//int mainCount(int argc, char* argv[]);
int mainCount(uint32_t numThreads, const std::string& indexBase, const LibraryFormat& libFmt,
              const std::vector<string>& undirReadFiles, const std::vector<string>& fwdReadFiles,
              const std::vector<string>& revReadFiles, const std::string& countFileOut,
              bool discardPolyA); 
int runIterativeOptimizer(int argc, char* argv[]);

int runKmerCounter(const std::string& sfCommand,
                   uint32_t numThreads,
                   const std::string& indexBase,
                   const LibraryFormat& libFmt,
                   const std::vector<string>& undirReadFiles,
                   const std::vector<string>& fwdReadFiles,
                   const std::vector<string>& revReadFiles,
                   const std::string& countFileOut,
                   bool discardPolyA) {

    /*
    std::stringstream argStream;

    argStream << sfCommand << " ";

    // if we're ignoring polyA/polyT in reads
    if (discardPolyA) { argStream << "--polya" << " "; }
    argStream << "--index " << indexBase << " ";
    argStream << "--counts " << countFileOut << " ";
    argStream << "--threads " << numThreads << " ";

    std::cerr << "Undirected reads:[";
    // Undirected reads
    if (undirReadFiles.size() > 0) {
        argStream << "--reads=";
        for (auto& rfile : undirReadFiles) {
            std::cerr << rfile << " ";
            argStream << rfile << " ";
        }
    }
    std::cerr << "]\n";

    std::cerr << "Sense reads:[";
    // Sense reads
    if (fwdReadFiles.size() > 0) {
      argStream << "--forward=";
      for (auto& rfile : fwdReadFiles) {
          std::cerr << rfile << " ";
          argStream << rfile << " ";
      }
    }
    std::cerr << "]\n";

    std::cerr << "Anti-sense reads:[";
    // Anti-sense reads
    if (revReadFiles.size() > 0) {
      argStream << "--reverse=";
      for (auto& rfile : revReadFiles) {
          std::cerr << rfile << " ";
          argStream << rfile << " ";
      }
    }
    std::cerr << "]\n";

    std::string argString = argStream.str();
    boost::trim(argString);
    */

    // Run the Sailfish counting phase as an separate process.
    // This will force the mmapped memory to be cleaned up.
    auto pid = fork();
    if (pid == 0) { //child
        /*
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
        */
        
        int ret = mainCount(numThreads, indexBase, libFmt, undirReadFiles,
                            fwdReadFiles, revReadFiles, countFileOut, discardPolyA); 
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
                          bool noBiasCorrect,
                          double minAbundance,
                          double maxDelta) {

  using std::vector;
  using std::string;
  using std::cerr;

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
    argStream << "--min_abundance " << minAbundance << " ";
    argStream << "--delta " << maxDelta << " ";
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

/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatString(std::string& fmt) {
    using std::vector;
    using std::string;
    using std::map;
    using std::stringstream;
 
    // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

    // first convert the string to upper-case
    for (auto& c : fmt) { c = std::toupper(c); }
    // split on the delimiter ':', and put the key, value (k=v) pairs into a map
    stringstream ss(fmt);
    string item;
    map<string, string> kvmap;
    while (std::getline(ss, item, ':')) {
        auto splitPos = item.find('=', 0);
        string key{item.substr(0, splitPos)};
        string value{item.substr(splitPos+1)};
        kvmap[key] = value;
    }
   
    map<string, ReadType> readType = {{"SE", ReadType::SINGLE_END}, {"PE", ReadType::PAIRED_END}};
    map<string, ReadOrientation> orientationType = {{">>", ReadOrientation::SAME},
                                           {"<>", ReadOrientation::AWAY},
                                           {"><", ReadOrientation::TOWARD},
                                           {"*", ReadOrientation::NONE}};
    map<string, ReadStrandedness> strandType = {{"SA", ReadStrandedness::SA}, 
                                    {"AS", ReadStrandedness::AS}, 
                                    {"A", ReadStrandedness::A}, 
                                    {"S", ReadStrandedness::S},
                                    {"U", ReadStrandedness::U}};
    auto it = kvmap.find("T");
    string typeStr = "";
    if (it != kvmap.end()) {
        typeStr = it->second;
    } else {
        it = kvmap.find("TYPE");
        if (it != kvmap.end()) {
            typeStr = it->second;
        }
    }
    
    if (typeStr != "SE" and typeStr != "PE") {
        string e = typeStr + " is not a valid read type; must be one of {SE, PE}";
        throw std::invalid_argument(e);
    }
   
    ReadType type = (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
    ReadOrientation orientation = (type == ReadType::SINGLE_END) ? ReadOrientation::NONE : ReadOrientation::TOWARD;
    ReadStrandedness strandedness{ReadStrandedness::U};
    // Construct the LibraryFormat class from the key, value map
    for (auto& kv : kvmap) {
        auto& k = kv.first; auto& v = kv.second;
        if (k == "O" or k == "ORIENTATION") {
            auto it = orientationType.find(v);
            if (it != orientationType.end()) { orientation = orientationType[it->first]; } else {
                string e = v + " is not a valid orientation type; must be one of {>>, <>, ><}";
                throw std::invalid_argument(e);
            }

        }
        if (k == "S" or k == "STRAND") {
            auto it = strandType.find(v);
            if (it != strandType.end()) { strandedness = strandType[it->first]; } else {
                string e = v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
                throw std::invalid_argument(e);
            }
        }
        
    }
    LibraryFormat lf(type, orientation, strandedness);
    return lf;
}

int mainQuantify( int argc, char *argv[] ) {

    using std::vector;
    using std::string;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    string sfCommand = argv[0];
    uint32_t maxThreads = std::thread::hardware_concurrency();
    bool noBiasCorrect = false;
    double minAbundance{0.01};
    double maxDelta{std::numeric_limits<double>::infinity()};

    std::vector<string> undirReadFiles;// = vm["reads"].as<std::vector<string>>();
    std::vector<string> fwdReadFiles;// = vm["forward"].as<std::vector<string>>();
    std::vector<string> revReadFiles;// = vm["reverse"].as<std::vector<string>>();

    std::vector<string> unmatedReadFiles;// = vm["reads"].as<std::vector<string>>();
    std::vector<string> mate1ReadFiles;// = vm["forward"].as<std::vector<string>>();
    std::vector<string> mate2ReadFiles;// = vm["reverse"].as<std::vector<string>>();



    // format string: (T|TYPE)=(SE|PE):(O|ORIENTATION)=(>>|<>|><):(S|STRAND)=(AS|SA|S|A|U)
    po::options_description generic("Sailfish quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "Sailfish index [output of the \"Sailfish index\" command")
    ("libtype,l", po::value<string>(), "Format string describing the library type")
    ("reads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
     "List of files containing reads of unknown orientation (i.e. \"undirected\" reads)")
    ("forward,F", po::value<vector<string>>(&fwdReadFiles)->multitoken(),
     "List of files containing reads oriented in the \"sense\" direction")
    ("reverse,R", po::value<vector<string>>(&revReadFiles)->multitoken(),
     "List of files containing reads oriented in the \"anti-sense\" direction")
    ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
        "File containing the #1 mates")
    ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
        "File containing the #2 mates")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
    ("min_abundance,m", po::value<double>(&minAbundance)->default_value(0.01),
     "transcripts with an abundance (KPKM) lower than this value will be reported at zero.")
    //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("out,o", po::value<string>(), "Basename of file where estimates are written")
    ("iterations,n", po::value<size_t>()->default_value(30), "number of iterations to run the optimzation")
    ("delta,d", po::value<double>(&maxDelta)->default_value(std::numeric_limits<double>::infinity()), "consider the optimization to have converged if the relative change in \n"
                                                           "the estimated abundance of all transcripts is below this threshold")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
    ("force,f", po::bool_switch(), "Force the counting phase to rerun, even if a count databse exists." )
    ("polya,a", po::bool_switch(), "polyA/polyT k-mers should be discarded")
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

        string libFmtStr = vm["libtype"].as<string>();
        LibraryFormat libFmt = parseLibraryFormatString(libFmtStr);

        if (libFmt.check()) {
            std::cerr << libFmt << "\n";
        } else {
            std::stringstream ss;
            ss << libFmt << " is invalid!";
            throw std::invalid_argument(ss.str());
        }
       
        // Ensure that we have only unmated reads with a single end library
        if (libFmt.type == ReadType::SINGLE_END) {
            if (unmatedReadFiles.size() == 0) {
                string e= "You must provide unmated read files with a single-end library type";
                throw std::invalid_argument(e);
            }
        }
        // or #1 and #2 mates with a paired-end library
        if (libFmt.type == ReadType::PAIRED_END) {
            if (mate1ReadFiles.size() == 0 or mate2ReadFiles.size() == 0) {
                string e = "You must provide #1 and #2 mated read files with a paired-end library type";
                throw std::invalid_argument(e);
            }
        }

        bfs::path indexBasePath(vm["index"].as<string>());
        bfs::path outputBasePath(vm["out"].as<string>());
        //string tgmap = vm["tgmap"].as<string>();
        //string tgmap = indexBase+".tgm";
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        bool force = vm["force"].as<bool>();
        bool discardPolyA = vm["polya"].as<bool>();

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
        
        // create the directory for log files
        bfs::path logDir = outputBasePath / "logs";
        boost::filesystem::create_directory(logDir);

        //string countFileOut = outputBase + ".sfc";
        bfs::path countFilePath(outputBasePath); countFilePath /= "reads.sfc";
        bfs::path indexPath(indexBasePath); indexPath /= "transcriptome";

        mustRecount = (force or !boost::filesystem::exists(countFilePath));
        if (mustRecount) {
            //          runKmerCounter(sfCommand, numThreads, indexPath.string(), undirReadFiles, fwdReadFiles, revReadFiles, countFilePath.string(), discardPolyA);
            runKmerCounter(sfCommand, numThreads, indexPath.string(), libFmt, unmatedReadFiles, mate1ReadFiles,
                           mate2ReadFiles, countFilePath.string(), discardPolyA);

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
                              iterations, lutBasePath, estFilePath,
                              noBiasCorrect, minAbundance, maxDelta);

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "============\n";
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << "============\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
    }

    return 0;
}
