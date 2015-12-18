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
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <chrono>
#include <iomanip>

#include "cereal/archives/binary.hpp"
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

//#include "KmerDist.hpp"
#include "SailfishUtils.hpp"
#include "SailfishConfig.hpp"
#include "VersionChecker.hpp"

int help(int argc, char* argv[]) {
  auto helpmsg = R"(
  ===============

  Please invoke sailfish with one of the following commands {index, quant, sf}.
  For more information on the options for theses particular methods, use the -h
  flag along with the method name.  For example:

  sailfish index -h

  will give you detailed help information about the index command.
  )";

  std::cerr << "  Sailfish v" << sailfish::version << helpmsg << "\n";
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

int mainIndex(int argc, char* argv[]);
int mainQuantify(int argc, char* argv[]);

bool verbose = false;

int main( int argc, char* argv[] ) {
  using std::string;
  namespace po = boost::program_options;

  try {

    po::options_description hidden("hidden");
    hidden.add_options()
    ("command", po::value<string>(), "command to run {index, quant, sf}");

    po::options_description sfopts("Allowed Options");
    sfopts.add_options()
    ("version,v", "print version string")
    ("no-version-check", "don't check with the server to see if this is the latest version")
    ("help,h", "produce help message")
    ;

    po::options_description all("Allowed Options");
    all.add(sfopts).add(hidden);

    // po::options_description sfopts("Command");
    // sfopts.add_options()

    po::positional_options_description pd;
    pd.add("command", 1);

    size_t topLevelArgc = argc;
    for (size_t i : boost::irange(size_t{1}, static_cast<size_t>(argc))) {
      if (argv[i][0] != '-') {
        topLevelArgc = i+1;
        break;
      }
    }

    po::variables_map vm;
    po::parsed_options parsed = po::command_line_parser(topLevelArgc, argv).options(all).positional(pd).allow_unregistered().run();
    po::store(parsed, vm);

/*    std::vector<string> subcommand_options = po::collect_unrecognized(parsed.options, po::include_positional);
    for (auto& s : subcommand_options) {
        std::cerr << "option: " << s << "\n";
    }
*/

    if (vm.count("version")) {
      std::cerr << "version : " << sailfish::version << "\n";
      std::exit(0);
    }

    if (vm.count("help") and !vm.count("command")) {
        std::cout << sfopts << std::endl;
        help(argc, argv);
        std::exit(0);
    }

    if (!vm.count("no-version-check")){
      std::string versionMessage = getVersionMessage();
      std::cerr << versionMessage;
    }

    po::notify(vm);

    std::unordered_map<string, std::function<int(int, char*[])>> cmds({
      {"index", mainIndex},
      {"quant", mainQuantify},
      {"sf", mainSailfish}
    });

    string cmd = vm["command"].as<string>();

    int subCommandArgc = argc - topLevelArgc + 1;
    char** argv2 = new char*[subCommandArgc];
    argv2[0] = argv[0];
    std::copy_n( &argv[topLevelArgc], argc-topLevelArgc, &argv2[1] );

    auto cmdMain = cmds.find(cmd);
    if (cmdMain == cmds.end()) {
      help(subCommandArgc, argv2);
    } else {
      cmdMain->second(subCommandArgc, argv2);
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
