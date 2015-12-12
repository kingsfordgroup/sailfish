#ifndef READ_KMER_DIST_HPP__
#define READ_KMER_DIST_HPP__
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
constexpr int64_t constExprPow(int64_t base, unsigned int exp)
{
	if(exp==0) return 1;
	else 
	{
		const int64_t res = constExprPow(base*base,exp/2);
		return (exp%2)?res*base:res;
	}
}

template <unsigned int K, typename CountT = uint32_t>
class ReadKmerDist {
public:
	std::array<CountT, constExprPow(4,K)> readCounts;
	ReadKmerDist():
	{
		for (size_t i = 0; i < readCounts.size(); ++i)
			readCounts[i] = 0; 
	}
	inline uint32_t indexForKmer(const char* s, bool reverseComplement) {
		// The index we'll return
		uint32_t idx{0};
		// The number of bits we need to shift the
		// current mask to the left.
		if(!reverseComplement)
		{
			for (size_t i = 0; i < K; ++i) {
				char c = s[i];
				switch (c) {
					case 'A':
					case 'a':
						break;
					case 'C':
					case 'c':
						idx += 1;
						break;
					case 'G':
					case 'g':
						idx += 2;
						break;
					case 'T':
					case 't':
					case 'U':
					case 'u':
						idx += 3;
						break;
					default:
						break;
				}
				idx = idx << 2;
			}
		}
		else
		{
			for(size_t i=K-1 ; i>=0 ; i--)
			{
				switch(s[i])
				{
					case 'T':
					case 't':
					case 'u':
					case 'U': break;
					case 'C':
					case 'c': idx += 2;
					break;
					case 'G':
					case 'g': idx += 1;
					break;
					case 'A':
					case 'a': idx += 3;
					break;
				}
				idx = idx << 2;
			}
		}
		return idx;
	}
	
	inline void updateReadCount(const char *s, const char *end,bool reverseComplement)
	{
		if (std::distance(s, end) >= K) {
			auto idx = indexForKmer(s,reverseComplement);
			readCounts[idx]++;
		}
	}
	
};
#endif 
