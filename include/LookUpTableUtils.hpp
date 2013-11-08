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


#ifndef LOOKUPTABLE_UTILS_HPP
#define LOOKUPTABLE_UTILS_HPP

#include <algorithm>
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

#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"

#include <boost/range/irange.hpp>
#include "ezETAProgressBar.hpp"

namespace LUTTools {

using TranscriptID = uint32_t;
using KmerID = uint64_t;
using ReadLength = int64_t;
using Offset = uint64_t;
using Length = uint32_t;
using TranscriptList = std::vector<TranscriptID>;

struct TranscriptInfo{
  TranscriptID transcriptID;
  TranscriptID geneID;
  std::string name;
  Length length;
  std::vector<KmerID> kmers; // TranscriptID => KmerID
};

/**
 *  \brief Dump the k-mer memberships vector to the file fname
 **/
void dumpKmerEquivClasses(
                          const std::vector<KmerID>& memberships,
                          const std::string& fname);

std::vector<KmerID> readKmerEquivClasses(const std::string& fname);

void dumpKmerLUT(
    std::vector<TranscriptList> &transcriptsForKmerClass,
    const std::string &fname);

void readKmerLUT(
    const std::string &fname,
    std::vector<TranscriptList> &transcriptsForKmer);


void writeTranscriptInfo (TranscriptInfo *ti, std::ofstream &ostream);

std::unique_ptr<TranscriptInfo> readTranscriptInfo(std::ifstream &istream);


std::vector<Offset> buildTLUTIndex(const std::string &tlutfname, size_t numTranscripts);


}

#endif //LOOKUPTABLE_UTILS_HPP
