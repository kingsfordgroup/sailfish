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


#include "PartitionRefiner.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>


PartitionRefiner::PartitionRefiner(LUTTools::KmerID numElem) :
numElem_(numElem),
membership_(std::vector<LUTTools::KmerID>(numElem+1, 0)),
maxSetIdx_(0)
{}


void PartitionRefiner::splitWith(std::vector<LUTTools::KmerID> splitter) {

  // If there is nothing in the splitter, then we shouldn't alter the equivalence classes;
  // just return
  if (splitter.size() == 0) { return; }

                // Sort the splitter so that we can build the # occurrence => elem map
                // efficiently
                std::sort(splitter.begin(), splitter.end());
                splitter.emplace_back(std::numeric_limits<size_t>::max());
                std::unordered_map<size_t, std::vector<size_t>> numOcc;
        /*
                std::cerr << "transcript: [";
                for (auto e : splitter){
                        std::cerr << e << ", ";
                }
                std::cerr << "]\n";
        */

                size_t i{1};
                size_t num{1};
                size_t e = splitter[0];
                size_t len = splitter.size();

                if (len == 1) { numOcc[1].emplace_back(e); }

                while ( i < len ) {
                        //std::cerr << "e = " << e << ", splitter[" << i << "] = " << splitter[i] << "\n";
                        if (splitter[i] == e) {
                                ++num;
                        } else {
                                numOcc[num].emplace_back(e);
                                num = 1;
                                e = splitter[i];
                        }
                ++i;
                }
        /*
                std::cerr << "{ ";
                for (auto& kv : numOcc) {
                        std::cerr << kv.first << " : [";
                        for (auto e : kv.second) {
                                std::cerr << e << ", ";
                        }
                        std::cerr << "], ";
                }
                std::cerr << "}\n";
        */
                struct IndexPair { size_t oldIndex; size_t newIndex; };

                for (auto& kv : numOcc) {
                        std::unordered_map<size_t, size_t> newSetIdx;
                        for (auto e : kv.second) {
                                auto currentSetIdx = membership_[e];
                                if (newSetIdx.find(currentSetIdx) == newSetIdx.end()) {
                                        newSetIdx[currentSetIdx] = ++maxSetIdx_;
                                }
                                membership_[e] = newSetIdx[membership_[e]];
                        }
                }

        }

void PartitionRefiner::relabel() {
  std::unordered_map<LUTTools::KmerID, LUTTools::KmerID> relabelMap;

  for (LUTTools::KmerID i = 0; i < membership_.size(); ++i) {
    auto e = membership_[i];
    if (relabelMap.find(e) == relabelMap.end()) {
      relabelMap[e] = static_cast<LUTTools::KmerID>(relabelMap.size());
    }
    membership_[i] = relabelMap[e];
  }

  std::cerr << "after relabling, there are " <<
                *std::max_element(membership_.begin(), membership_.end())+1 << " eq classes\n";
}

const std::vector<LUTTools::KmerID>& PartitionRefiner::partitionMembership() { return membership_; }
