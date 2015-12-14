#ifndef READ_KMER_DIST_HPP__
#define READ_KMER_DIST_HPP__
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include "UtilityFunctions.hpp"

template <uint32_t K, typename CountT = uint32_t>
class ReadKmerDist {
public:
	std::array<CountT, constExprPow(4,K)> counts;

	ReadKmerDist() {
        // set a pseudo-count of 1
		for (size_t i = 0; i < counts.size(); ++i) {
			counts[i] = 1;
        }
	}

    inline uint32_t getK() { return K; }

    inline uint64_t totalCount() {
        CountT c{0};
        for (auto const& rc : counts) { c += rc; }
        return c;
    }

    // update the k-mer context for the hit at position p.
    // The underlying transcript is from [start, end)
	inline bool update(const char* start, const char *p, const char *end, bool isForward) {
        int posBeforeHit = 2;
        int posAfterHit = 4;
        bool success{false};
        if (isForward) {
            // If we can fit the window before and after the read
            if ((p - start) >= posBeforeHit and
                ((p - posBeforeHit + K) < end) ) {
                p -= posBeforeHit;
                // If the read matches in the forward direction, we take
                // the RC sequence.
                auto idx = indexForKmer(p, K, true);
                if (idx > counts.size()) {
                    std::cerr << "start! = " << *start << '\n';
                    std::cerr << "p = " << *p << '\n';
                    std::cerr << "p+1 = " << *(p+1) << '\n';
                    std::cerr << "p+2 = " << *(p+2) << '\n';
                    std::cerr << "p+3 = " << *(p+3) << '\n';
                    std::cerr << "p+4 = " << *(p+4) << '\n';
                    std::cerr << "p+5 = " << *(p+5) << '\n';
                    std::cerr << "idx = " << idx << ", size = " << counts.size() << '\n';
                }
                counts[idx]++;
                success = true;
            }
        } else {
            if ((p - start) >= posAfterHit and
                ((p - posAfterHit + K) < end) ) {
                p -= posAfterHit;
                auto idx = indexForKmer(p, K, false);
                if (idx > counts.size()) {
                    std::cerr << "start! = " << *start << '\n';
                    std::cerr << "p = " << *p << '\n';
                    std::cerr << "p+1 = " << *(p+1) << '\n';
                    std::cerr << "p+2 = " << *(p+2) << '\n';
                    std::cerr << "p+3 = " << *(p+3) << '\n';
                    std::cerr << "p+4 = " << *(p+4) << '\n';
                    std::cerr << "p+5 = " << *(p+5) << '\n';
                    std::cerr << "idx = " << idx << ", size = " << counts.size() << '\n';
                }
                counts[idx]++;
                success = true;
            }
		}
        return success;
	}

};
#endif
