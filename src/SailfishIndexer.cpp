/**
>HEADER
    Copyright (c) 2015 Rob Patro rob.patro@cs.stonybrook.edu

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

#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <cassert>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "jellyfish/config.h"
#include "jellyfish/err.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/jellyfish.hpp"

#include "cereal/archives/binary.hpp"

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#include "RapMapSAIndex.hpp"
#include "SailfishUtils.hpp"
#include "SailfishIndex.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

int mainIndex( int argc, char *argv[] ) {
    using std::string;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();
    bool useStreamingParser = false;

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<std::vector<string>>()->multitoken()->required(), "Transcript fasta file(s)." )
    //("tgmap,m", po::value<string>(), "file that maps transcripts to genes")
    ("kmerSize,k", po::value<uint32_t>()->required()->default_value(31), "Kmer size.")
    ("out,o", po::value<string>()->required(), "Output stem [all files needed by Sailfish will be of the form stem.*].")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use concurrently.")
    ("force,f", po::bool_switch(), "" )
    ;

    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
index
==========
Builds a Sailfish index
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        uint32_t merLen = vm["kmerSize"].as<uint32_t>();
        string outputStem = vm["out"].as<string>();
        std::vector<string> transcriptFiles = vm["transcripts"].as<std::vector<string>>();
        uint32_t numThreads = vm["threads"].as<uint32_t>();
        bool force = vm["force"].as<bool>();

        // Check to make sure that the specified output directory either doesn't exist, or is
        // a valid path (e.g. not a file)
        namespace bfs = boost::filesystem;

        // Ensure that the transcript file provided by the user exists
        auto txpFile = transcriptFiles.front();
        if (!bfs::exists(txpFile)) {
            std::cerr << "The provided transcript file [" << txpFile << "] does not seem to exist!\n";
            std::cerr << "Please check that the correct path was provided.\n";
            std::exit(1);
        }
        // and that it is, in fact, a file
        if (bfs::is_directory(txpFile)) {
            std::cerr << "The provided transcript file [" << txpFile << "] appears to be a directory!\n";
            std::cerr << "Please check that the correct path was provided.\n";
            std::exit(1);
        }

        // Check that the output path doesn't exist yet (or at least is not a file)
        if (bfs::exists(outputStem) and !bfs::is_directory(outputStem)) {
            std::cerr << "The provided output path [" << outputStem << "] " <<
                         "already exists and is not a directory\n.";
            std::cerr << "Please either provide a different path or " <<
                         "delete the existing file.\n";
            std::exit(1);
        }

        bool mustRecompute = false;
        if (!bfs::exists(outputStem)) {
            mustRecompute = true;
            try {
                bool success = bfs::create_directories(outputStem);
                if (!success) { throw std::runtime_error("unspecified error creating file."); }
            } catch ( std::exception& e ) {
                std::cerr << "Creation of " << outputStem << " failed [" << e.what() << "]\n.";
                std::cerr << "Exiting.\n";
                return 1;
            }
        }

        bfs::path outputPath(outputStem);

        // create the directory for log files
        bfs::path logDir = outputPath / "logs";
        boost::filesystem::create_directories(logDir);

        bfs::path logPath = logDir / "sailfish_index.log";
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        //auto consoleLog = spdlog::create("consoleLog", {consoleSink});
        //auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        std::cerr << "writing log to " << logPath.string() << "\n";

        /** No more of this bias correction method!
        // First, compute the transcript features in case the user
        // ever wants to bias-correct his / her results
        bfs::path transcriptBiasFile(outputPath); transcriptBiasFile /= "bias_feats.txt";

        jointLog->info() << "computeBiasFeatures( {";
        for (auto& tf : transcriptFiles) {
            jointLog->info() << "[" << tf << "] ";
        }
   	    std::string txpBiasFileName = transcriptBiasFile.string();
        jointLog->info(", {}, {}, {}\n", txpBiasFileName, useStreamingParser, numThreads);
        computeBiasFeatures(transcriptFiles, transcriptBiasFile, useStreamingParser, numThreads);
        */

        bfs::path headerPath = outputPath / "header.json";
        mustRecompute = (force or !boost::filesystem::exists(headerPath));

        if (!mustRecompute) {
            // Check that the jellyfish has at the given location
            // was computed with the correct kmer length.
            jointLog->info("Index exists but will not be rebuilt --- use the force "
                           "option to rebuild the index");
        }

        if (mustRecompute) {

            std::vector<const char*> argVec;
            argVec.push_back("foo");
            argVec.push_back("-k");

            if (merLen % 2 == 0) {
                jointLog->error("k-mer length should be odd to avoid a k-mer being it's own reverse complement\n"
                                "please specify an odd value of k\n");
                jointLog->flush();
                spdlog::drop("jointLog");
                return 1;
            }

            fmt::MemoryWriter optWriter;
            optWriter << merLen;
            argVec.push_back(optWriter.str().c_str());
            argVec.push_back("-t");
            argVec.push_back(transcriptFiles.front().c_str());
            argVec.push_back("-i");
            argVec.push_back(outputPath.string().c_str());
            SailfishIndex sidx(jointLog);
            sidx.build(outputPath, argVec, merLen);
        } else {
            std::cerr << "All index files seem up-to-date.\n";
            std::cerr << "To force Sailfish to rebuild the index, use the --force option.\n";
        }

    } catch (po::error &e) {
        std::cerr << "Exception: [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception: [" << e.what() << "]\n";
        std::cerr << argv[0] << " index was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }

    return 0;
}
