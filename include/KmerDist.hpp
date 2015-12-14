/**
 > HEADER *               *
 Copyright (c) 2015 Rob Patro rob.patro@cs.stonybrook.edu
 This file is part of Sailfish / Salmon.
 Sailfish is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 Sailfish is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 You should have received a copy of the GNU General Public License
 along with Sailfish. If not, see <http://www.gnu.org/licenses/>.
 <HEADER
 **/
#ifndef KMER_DIST_HPP__
#define KMER_DIST_HPP__
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdint>
#include "UtilityFunctions.hpp"

// templatized on the k-mer length
template <unsigned int K, typename CountT=uint32_t>
class KmerDist {
	using Bias = double;
	using Offset = size_t;
public:
	// The number of k-mers we'll have to store given the k-mer size
	//constexpr const uint64_t numKmers{constExprPow(4, K)};
	std::array<CountT, constExprPow(4,K)> counts;
	std::map<uint64_t,bool> hexamers;

	const uint64_t numPossKmers = constExprPow(4,K);

	KmerDist() : haveBias_(false) {
		// set a pseudo-count for each k-mer
		for (size_t i = 0; i < counts.size(); ++i) { counts[i] = 0; }
	}


	inline void addKmers(const char* s, const char* end, bool reverseComplement) {
		// Make sure there are at least k bases left

/*		auto idx = indexForKmer(start,reverseComplement);
		char *it = start+1;
		while(std::distance(start, end) >= K) {
			auto idx = indexForKmer(start,reverseComplement);
			auto it = hexamers.find(idx);
			if(it==hexamers.end())
				hexamers[idx] = true;
			counts[idx]++;
		}
		return true;
*/
		if(!reverseComplement) {
			auto idx = indexForKmer(start, K, false);
			hexamers[idx]=true;
			counts[idx]++;
			uint32_t i=1;
			while(std::distance(start+i,end)>=K) {

				idx = nextKmerIndex(idx, start+i, K, false);
				auto it = hexamers.find(idx);
				if(it==hexamers.end())
					hexamers[idx] = true;
				counts[idx]++;
				i++;

			}
		} else {
			auto idx = indexForKmer(end-1, K, true);
			hexamers[idx]=true;
			counts[idx]++;
			uint32_t i=2;
			while(end-i>=start)
			{
				idx = nextKmerIndex(idx, end-i, K, true);
				auto it = hexamers.find(idx);
				if(it==hexamers.end())
					hexamers[idx] = true;
				counts[idx]++;
				i--;

			}
		}
	}
private:
	bool haveBias_;
	//std::vector<Bias> biases_;
	//std::unordered_map<std::string, size_t> offsetMap_;
};
#endif // KMER_DIST_HPP__
