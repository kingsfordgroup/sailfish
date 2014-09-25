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


#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <sstream>
#include <exception>
#include <random>
#include <queue>
#include <unordered_map>
#include "btree_map.h"
#include "btree_set.h"


/** BWA Includes */
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cctype>

// C++ string formatting library
#include "format.h"

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
}

#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

/** Boost Includes */
#include <boost/filesystem.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/partitioner.h"

#if HAVE_LOGGER
#include "g2logworker.h"
#include "g2log.h"
#endif

#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include "ClusterForest.hpp"
#include "PerfectHashIndex.hpp"
#include "LookUpTableUtils.hpp"
#include "SailfishMath.hpp"
#include "Transcript.hpp"
#include "LibraryFormat.hpp"
#include "SailfishUtils.hpp"
#include "SalmonUtils.hpp"
#include "ReadLibrary.hpp"
#include "SalmonConfig.hpp"

#include "PairSequenceParser.hpp"
#include "FragmentList.hpp"
#include "FragmentLengthDistribution.hpp"
#include "ReadExperiment.hpp"

extern unsigned char nst_nt4_table[256];
char* bwa_pg = "cha";

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;

using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

uint32_t transcript(uint64_t enc) {
    uint32_t t = (enc & 0xFFFFFFFF00000000) >> 32;
    return t;
}

uint32_t offset(uint64_t enc) {
    uint32_t o = enc & 0xFFFFFFFF;
    return o;
}

class Alignment {
    public:
        Alignment(TranscriptID transcriptIDIn, uint32_t kCountIn = 1, uint32_t fragLengthIn= 0, double logProbIn = sailfish::math::LOG_0) :
            transcriptID_(transcriptIDIn), kmerCount(kCountIn), fragLength_(fragLengthIn), logProb(logProbIn) {}

        inline TranscriptID transcriptID() { return transcriptID_; }
        inline uint32_t fragLength() { return fragLength_; }

        uint32_t kmerCount;
        double logProb;

    private:
        TranscriptID transcriptID_;
        uint32_t fragLength_;
};

void processMiniBatch(
        double logForgettingMass,
        std::vector<std::vector<Alignment>>& batchHits,
        std::vector<Transcript>& transcripts,
        ClusterForest& clusterForest,
        FragmentLengthDistribution& fragLengthDist,
        std::atomic<uint64_t>& numAssignedFragments,
        std::default_random_engine& randEng,
        bool initialRound,
        bool& burnedIn
        ) {

    using sailfish::math::LOG_0;
    using sailfish::math::LOG_1;
    using sailfish::math::logAdd;
    using sailfish::math::logSub;

    size_t numTranscripts{transcripts.size()};
    size_t numReads{batchHits.size()};
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

    bool updateCounts = initialRound;

    // Build reverse map from transcriptID => hit id
    using HitID = uint32_t;
    btree::btree_map<TranscriptID, std::vector<Alignment*>> hitsForTranscript;
    size_t hitID{0};
    for (auto& hv : batchHits) {
        for (auto& tid : hv) {
            hitsForTranscript[tid.transcriptID()].push_back(&tid);
        }
        ++hitID;
    }

    double clustTotal = std::log(batchHits.size()) + logForgettingMass;
    {
        // E-step

        // Iterate over each group of alignments (a group consists of all alignments reported
        // for a single read).  Distribute the read's mass proportionally dependent on the
        // current hits
        for (auto& alnGroup : batchHits) {
            if (alnGroup.size() == 0) { continue; }
            double sumOfAlignProbs{LOG_0};
            // update the cluster-level properties
            bool transcriptUnique{true};
            auto firstTranscriptID = alnGroup.front().transcriptID();
            std::unordered_set<size_t> observedTranscripts;
            for (auto& aln : alnGroup) {
                auto transcriptID = aln.transcriptID();
                auto& transcript = transcripts[transcriptID];
                transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;

                double logFragProb = sailfish::math::LOG_1;

                // The alignment probability is the product of a transcript-level term (based on abundance and) an alignment-level
                // term below which is P(Q_1) * P(Q_2) * P(F | T)
                double logRefLength = std::log(refLength);
                double transcriptLogCount = transcript.mass();
                if ( transcriptLogCount != LOG_0 ) {
                    double errLike = sailfish::math::LOG_1;
                    if (burnedIn) {
                        //errLike = errMod.logLikelihood(aln, transcript);
                    }

                    double logFragProb = (aln.fragLength() == 0) ? LOG_1 : fragLengthDist.pmf(static_cast<size_t>(aln.fragLength()));

                    //aln.logProb = std::log(aln.kmerCount) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    //aln.logProb = std::log(std::pow(aln.kmerCount,2.0)) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    aln.logProb = (transcriptLogCount - logRefLength) + logFragProb;// + qualProb + errLike;

                    sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);

                    if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                        if (updateCounts) { transcripts[transcriptID].addTotalCount(1); }
                        observedTranscripts.insert(transcriptID);
                    }
                } else {
                    aln.logProb = LOG_0;
                }
            }

            // normalize the hits
            if (sumOfAlignProbs == LOG_0) { std::cerr << "0 probability fragment; skipping\n"; continue; }
            for (auto& aln : alnGroup) {
                aln.logProb -= sumOfAlignProbs;
                auto transcriptID = aln.transcriptID();
                auto& transcript = transcripts[transcriptID];

                double r = uni(randEng);
                if (!burnedIn and r < std::exp(aln.logProb)) {
                    //errMod.update(aln, transcript, aln.logProb, logForgettingMass);
                    double fragLength = aln.fragLength();
                    if (fragLength > 0.0) {
                        //if (aln.fragType() == ReadType::PAIRED_END) {
                        fragLengthDist.addVal(fragLength, logForgettingMass);
                    }
                }


            } // end normalize

            // update the single target transcript
            if (transcriptUnique) {
                if (updateCounts) {
                    transcripts[firstTranscriptID].addUniqueCount(1);
                }
                clusterForest.updateCluster(firstTranscriptID, 1, logForgettingMass, updateCounts);
            } else { // or the appropriate clusters
                clusterForest.mergeClusters<Alignment>(alnGroup.begin(), alnGroup.end());
                clusterForest.updateCluster(alnGroup.front().transcriptID(), 1, logForgettingMass, updateCounts);
            }

            } // end read group
        }// end timer

        double individualTotal = LOG_0;
        {
            // M-step
            double totalMass{0.0};
            for (auto kv = hitsForTranscript.begin(); kv != hitsForTranscript.end(); ++kv) {
                auto transcriptID = kv->first;
                // The target must be a valid transcript
                if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr << "index " << transcriptID << " out of bounds\n"; }

                auto& transcript = transcripts[transcriptID];

                // The prior probability
                double hitMass{LOG_0};

                // The set of alignments that match transcriptID
                auto& hits = kv->second;
                std::for_each(hits.begin(), hits.end(), [&](Alignment* aln) -> void {
                        if (!std::isfinite(aln->logProb)) { std::cerr << "hitMass = " << aln->logProb << "\n"; }
                        hitMass = logAdd(hitMass, aln->logProb);
                        });

                double updateMass = logForgettingMass + hitMass;
                individualTotal = logAdd(individualTotal, updateMass);
                totalMass = logAdd(totalMass, updateMass);
                transcript.addMass(updateMass);
            } // end for
        } // end timer
        numAssignedFragments += numReads;
        if (numAssignedFragments >= 5000000 and !burnedIn) { burnedIn = true; }
}

uint32_t basesCovered(std::vector<uint32_t>& kmerHits) {
    std::sort(kmerHits.begin(), kmerHits.end());
    uint32_t covered{0};
    uint32_t lastHit{0};
    uint32_t kl{20};
    for (auto h : kmerHits) {
        covered += std::min(h - lastHit, kl);
        lastHit = h;
    }
    return covered;
}

uint32_t basesCovered(std::vector<uint32_t>& posLeft, std::vector<uint32_t>& posRight) {
    return basesCovered(posLeft) + basesCovered(posRight);
}

class KmerVote {
    public:
        KmerVote(uint32_t vp, uint32_t rp, uint32_t vl) : votePos(vp), readPos(rp), voteLen(vl) {}
        uint32_t votePos{0};
        uint32_t readPos{0};
        uint32_t voteLen{0};
};

class MatchFragment {
    public:
        MatchFragment(uint32_t refStart_, uint32_t queryStart_, uint32_t length_) :
            refStart(refStart_), queryStart(queryStart_), length(length_) {}

        uint32_t refStart, queryStart, length;
        uint32_t weight;
        double score;
};

bool precedes(const MatchFragment& a, const MatchFragment& b) {
    return (a.refStart + a.length) < b.refStart and
           (a.queryStart + a.length) < b.queryStart;
}


class TranscriptHitList {
    public:
        uint32_t bestHitPos{0};
        uint32_t bestHitScore{0};

        std::vector<KmerVote> votes;
        std::vector<KmerVote> rcVotes;

        uint32_t targetID;


        void addFragMatch(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            uint32_t votePos = (readPos > tpos) ? 0 : tpos - readPos;
            votes.emplace_back(votePos, readPos, voteLen);
        }

        void addFragMatchRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            uint32_t votePos = (readPos > tpos) ? 0 : tpos - readPos;
            rcVotes.emplace_back(votePos, readPos, voteLen);
        }

        uint32_t totalNumHits() { return std::max(votes.size(), rcVotes.size()); }

        void computeBestLoc_(std::vector<KmerVote>& sVotes, uint32_t& maxClusterPos, uint32_t& maxClusterCount) {
            if (sVotes.size() == 0) { return; }
            struct VoteInfo {
                uint32_t coverage;
                uint32_t rightmostBase;
            };

            boost::container::flat_map<uint32_t, VoteInfo> hitMap;
            int32_t currClust{static_cast<int32_t>(sVotes.front().votePos)};
            for (size_t j = 0; j < sVotes.size(); ++j) {

                uint32_t votePos = sVotes[j].votePos;
                uint32_t readPos = sVotes[j].readPos;
                uint32_t voteLen = sVotes[j].voteLen;

                if (votePos >= currClust) {
                    if (votePos - currClust > 10) {
                        currClust = votePos;
                    }
                    auto& hmEntry = hitMap[currClust];

                    hmEntry.coverage += std::min(voteLen, (votePos + readPos + voteLen) - hmEntry.rightmostBase);
                    hmEntry.rightmostBase = votePos + readPos + voteLen;
                } else if (votePos < currClust) {
                    std::cerr << "Should not have votePos = " << votePos << " <  currClust = " << currClust << "\n";
                    std::exit(1);
                }

                if (hitMap[currClust].coverage > maxClusterCount) {
                    maxClusterCount = hitMap[currClust].coverage;
                    maxClusterPos = currClust;
                }

            }

        }

        bool computeBestChain() {
            std::sort(votes.begin(), votes.end(),
                    [](const KmerVote& v1, const KmerVote& v2) -> bool {
                        if (v1.votePos == v2.votePos) {
                            return v1.readPos < v2.readPos;
                        }
                        return v1.votePos < v2.votePos;
                    });

            std::sort(rcVotes.begin(), rcVotes.end(),
                    [](const KmerVote& v1, const KmerVote& v2) -> bool {
                        if (v1.votePos == v2.votePos) {
                            return v1.readPos < v2.readPos;
                        }
                        return v1.votePos < v2.votePos;
                    });

            uint32_t maxClusterPos{0};
            uint32_t maxClusterCount{0};

            computeBestLoc_(votes, maxClusterPos, maxClusterCount);
            computeBestLoc_(rcVotes, maxClusterPos, maxClusterCount);

            bestHitPos = maxClusterPos;
            bestHitScore = maxClusterCount;

            return true;
        }

            /*
        std::vector<MatchFragment> frags;
        std::vector<MatchFragment> rcFrags;

        void addFragMatch(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            frags.emplace_back(tpos, readPos, voteLen);
        }

        void addFragMatchRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
             rcFrags.emplace_back(tpos, readPos, voteLen);
        }

        double gapCost_(MatchFragment& a, MatchFragment& b) {
            return 0.5 * ((a.queryStart + a.length - b.queryStart) +
                          (a.refStart + a.length - b.refStart) );
        }

        void computeBestChain_(std::vector<MatchFragment>& mfrags, uint32_t& maxClusterPos, uint32_t& maxClusterCount) {
            if (mfrags.size() == 0) { return; }
            std::vector<double> scores(mfrags.size(), 0);

            for (size_t i = 0; i < mfrags.size(); ++i) {
                auto& fragI = mfrags[i];
                for (size_t j = 0; j < mfrags.size(); ++j) {
                    auto& fragJ = mfrags[j];
                    if (precedes(fragJ, fragI)) {
                        auto score = scores[j] - gapCost_(fragJ, fragI);
                        if (score > scores[i]) { scores[i] = score; }
                    }
                }
                scores[i] += fragI.weight;
                if (scores[i] > maxClusterCount) {
                    maxClusterCount = scores[i];
                    maxClusterPos = i;
                }
            }
        }

        bool computeBestChain() {
            uint32_t maxClusterPos{0};
            uint32_t maxClusterCount{0};

            computeBestChain_(frags, maxClusterPos, maxClusterCount);
            computeBestChain_(rcFrags, maxClusterPos, maxClusterCount);

            bestHitPos = maxClusterPos;
            bestHitScore = maxClusterCount;

            return true;
        }
        */
};

template <typename CoverageCalculator>
void getHitsForFragment(std::pair<header_sequence_qual, header_sequence_qual>& frag,
                        bwaidx_t *idx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        int minLen,
                        int minIWidth,
                        int splitWidth,
                        double coverageThresh,
                        std::vector<Alignment>& hitList,
                        uint64_t& hitListCount,
                        std::vector<Transcript>& transcripts) {

    uint64_t leftHitCount{0};

    //std::unordered_map<uint64_t, TranscriptHitList> leftHits;
    //std::unordered_map<uint64_t, TranscriptHitList> rightHits;

    //std::unordered_map<uint64_t, CoverageCalculator> leftHits;
    //std::unordered_map<uint64_t, CoverageCalculator> rightHits;


    std::unordered_map<uint64_t, CoverageCalculator> leftHits;
    std::unordered_map<uint64_t, CoverageCalculator> rightHits;
    //leftHits.set_empty_key(std::numeric_limits<uint64_t>::max());
    //rightHits.set_empty_key(std::numeric_limits<uint64_t>::max());

    uint32_t leftReadLength{0};
    uint32_t rightReadLength{0};

    //---------- Left End ----------------------//
    {
        std::string& readStr   = frag.first.seq;
        uint32_t readLen      = frag.first.seq.size();

        leftReadLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        char* readPtr = const_cast<char*>(readStr.c_str());

        smem_set_query(itr, readLen, reinterpret_cast<uint8_t*>(readPtr));

        // while there are more matches on the query
        while ((a = smem_next(itr)) != 0) {
        //while ((a = smem_next(itr, minLen<<1, splitWidth)) != 0) {

            for (size_t mi = 0; mi < a->n; ++mi) {
                bwtintv_t *p = &a->a[mi];
                if (static_cast<uint32_t>(p->info) - (p->info>>32) < minLen) continue;
                uint32_t qstart = static_cast<uint32_t>(p->info>>32);
                uint32_t qend = static_cast<uint32_t>(p->info);
                // long numHits = static_cast<long>(p->x[2]);

                if ( p->x[2] <= minIWidth) {
                    for (bwtint_t k = 0; k < p->x[2]; ++k) {
                        bwtint_t pos;
                        int len, isRev, refID;
                        len = static_cast<uint32_t>(p->info - (p->info >> 32));
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &isRev);
                        if (isRev) { pos -= len - 1; }
                        bns_cnt_ambi(idx->bns, pos, len, &refID);
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset) + 1;

                        if (isRev) {
                            qstart = leftReadLength - qend;
                            leftHits[refID].addFragMatchRC(hitLoc, qstart, len);
                        } else {
                            leftHits[refID].addFragMatch(hitLoc, qstart, len);
                        }
                    } // for k
                } // if <= minIWidth

            } // for mi
        } // for all query matches
    }

    //---------- Right End ----------------------//
    {
        std::string& readStr   = frag.second.seq;
        uint32_t readLen      = frag.second.seq.size();

        rightReadLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        char* readPtr = const_cast<char*>(readStr.c_str());
        smem_set_query(itr, readLen, reinterpret_cast<uint8_t*>(readPtr));

        // while there are more matches on the query
        while ((a = smem_next(itr)) != 0) {
        //while ((a = smem_next(itr, minLen<<1, splitWidth)) != 0) {
            for (size_t mi = 0; mi < a->n; ++mi) {
                bwtintv_t *p = &a->a[mi];
                if (static_cast<uint32_t>(p->info) - (p->info>>32) < minLen) continue;
                uint32_t qstart = static_cast<uint32_t>(p->info>>32);
                uint32_t qend = static_cast<uint32_t>(p->info);
                // long numHits = static_cast<long>(p->x[2]);

                if ( p->x[2] <= minIWidth) {
                    for (bwtint_t k = 0; k < p->x[2]; ++k) {
                        bwtint_t pos;
                        int len, isRev, refID;
                        len = static_cast<uint32_t>(p->info - (p->info >> 32));
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &isRev);
                        if (isRev) { pos -= len - 1; }
                        bns_cnt_ambi(idx->bns, pos, len, &refID);
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset) + 1;

                        if (isRev) {
                            qstart = rightReadLength - qend;
                            rightHits[refID].addFragMatchRC(hitLoc, qstart, len);
                        } else {
                            rightHits[refID].addFragMatch(hitLoc, qstart, len);
                        }
                    } // for k
                } // if <= minIWidth

            } // for mi
        } // for all query matches

    } // end right

    size_t readHits{0};
    hitList.clear();

    for (auto& tHitList : leftHits) {
        // Coverage score
        tHitList.second.computeBestChain();
        ++leftHitCount;
    }

    double cutoffLeft{ coverageThresh * leftReadLength};
    double cutoffRight{ coverageThresh * rightReadLength};

    for (auto& tHitList : rightHits) {
        auto it = leftHits.find(tHitList.first);
        // Coverage score
        if (it != leftHits.end() and it->second.bestHitScore >= cutoffLeft) {
            tHitList.second.computeBestChain();
            if (tHitList.second.bestHitScore < cutoffRight) { continue; }
            uint32_t score = it->second.bestHitScore + tHitList.second.bestHitScore;
            uint32_t fragLength = std::abs(static_cast<int32_t>(it->second.bestHitPos) -
                                           static_cast<int32_t>(tHitList.second.bestHitPos +
                                                                rightReadLength));
            hitList.emplace_back(tHitList.first, score, fragLength);
            readHits += score;
            ++hitListCount;
        }
    }

}

/**
  *   Get hits for single-end fragment
  *
  *
  */
template <typename CoverageCalculator>
void getHitsForFragment(jellyfish::header_sequence_qual& frag,
                        bwaidx_t *idx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        int minLen,
                        int minIWidth,
                        int splitWidth,
                        double coverageThresh,
                        std::vector<Alignment>& hitList,
                        uint64_t& hitListCount,
                        std::vector<Transcript>& transcripts) {

    uint64_t leftHitCount{0};

    //std::unordered_map<uint64_t, TranscriptHitList> hits;
    std::unordered_map<uint64_t, CoverageCalculator> hits;

    uint32_t readLength{0};

    //---------- get hits ----------------------//
    {
        std::string& readStr   = frag.seq;
        uint32_t readLen      = frag.seq.size();

        readLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        char* readPtr = const_cast<char*>(readStr.c_str());

        smem_set_query(itr, readLen, reinterpret_cast<uint8_t*>(readPtr));

        // while there are more matches on the query
        while ((a = smem_next(itr)) != 0) {
        //while ((a = smem_next(itr, minLen<<1, splitWidth)) != 0) {
            for (size_t mi = 0; mi < a->n; ++mi) {
                bwtintv_t *p = &a->a[mi];
                if (static_cast<uint32_t>(p->info) - (p->info>>32) < minLen) continue;
                uint32_t qstart = static_cast<uint32_t>(p->info>>32);
                uint32_t qend = static_cast<uint32_t>(p->info);
                // long numHits = static_cast<long>(p->x[2]);

                if ( p->x[2] <= minIWidth) {
                    for (bwtint_t k = 0; k < p->x[2]; ++k) {
                        bwtint_t pos;
                        int len, isRev, refID;
                        len = static_cast<uint32_t>(p->info - (p->info >> 32));
                        pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &isRev);
                        if (isRev) { pos -= len - 1; }
                        bns_cnt_ambi(idx->bns, pos, len, &refID);
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset) + 1;
                        if (isRev) {
                            qstart = readLength - qend;
                            hits[refID].addFragMatchRC(hitLoc, qstart, len);
                        } else {
                            hits[refID].addFragMatch(hitLoc, qstart, len);
                        }
                    } // for k
                } // if <= minIWidth

            } // for mi
        } // for all query matches
    }

    size_t readHits{0};
    hitList.clear();

    double cutoff{ coverageThresh * readLength};
    for (auto& tHitList : hits) {
        // Coverage score
        tHitList.second.computeBestChain();
        if (tHitList.second.bestHitScore >= cutoff) {
            uint32_t score = tHitList.second.bestHitScore;
            hitList.emplace_back(tHitList.first, score);
            readHits += score;
            ++hitListCount;
            ++leftHitCount;
        }
    }

}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT, typename CoverageCalculator>
void processReadsMEM(ParserT* parser,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               bwaidx_t *idx,
               std::vector<Transcript>& transcripts,
               std::atomic<uint64_t>& batchNum,
               double& logForgettingMass,
               std::mutex& ffMutex,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               uint32_t minMEMLength,
               uint32_t maxMEMOcc,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn
               ) {
  uint64_t count_fwd = 0, count_bwd = 0;

  double forgettingFactor{0.65};

  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  std::vector<std::vector<Alignment>> hitLists;
  hitLists.resize(5000);

  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  // Super-MEM iterator
  smem_i *itr = smem_itr_init(idx->bwt);
  const bwtintv_v *a = nullptr;
  int minLen{static_cast<int>(minMEMLength)};
  int minIWidth{static_cast<int>(maxMEMOcc)};
  int splitWidth{0};

  size_t locRead{0};
  while(true) {
    typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    hitLists.resize(j->nb_filled);
    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch

        getHitsForFragment<CoverageCalculator>(j->data[i], idx, itr, a,
                                               minLen, minIWidth, splitWidth,
                                               coverageThresh,
                                               hitLists[i], hitListCount,
                                               transcripts);
        auto& hitList = hitLists[i];

        locRead++;
        ++numObservedFragments;
        if (numObservedFragments % 50000 == 0) {
    	    iomutex.lock();
            const char RESET_COLOR[] = "\x1b[0m";
            char green[] = "\x1b[30m";
            green[3] = '0' + static_cast<char>(fmt::GREEN);
            char red[] = "\x1b[30m";
            red[3] = '0' + static_cast<char>(fmt::RED);
            fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red, numObservedFragments, green, RESET_COLOR);
    	    iomutex.unlock();
        }

        // If the read mapped to > 100 places, discard it
        if (hitList.size() > 100) { hitList.clear(); }

    } // end for i < j->nb_filled

    auto oldBatchNum = batchNum++;
    if (oldBatchNum > 1) {
        ffMutex.lock();
        logForgettingMass += forgettingFactor * std::log(static_cast<double>(oldBatchNum-1)) -
        std::log(std::pow(static_cast<double>(oldBatchNum), forgettingFactor) - 1);
        ffMutex.unlock();
    }

    processMiniBatch(logForgettingMass, hitLists, transcripts, clusterForest,
                     fragLengthDist, numAssignedFragments, eng, initialRound, burnedIn);
  }
  smem_itr_destroy(itr);
}

int performBiasCorrection(boost::filesystem::path featPath,
                          boost::filesystem::path expPath,
                          double estimatedReadLength,
                          double kmersPerRead,
                          uint64_t mappedKmers,
                          uint32_t merLen,
                          boost::filesystem::path outPath,
                          size_t numThreads);

void processReadLibrary(
        ReadLibrary& rl,
        bwaidx_t* idx,
        std::vector<Transcript>& transcripts,
        ClusterForest& clusterForest,
        std::atomic<uint64_t>& numObservedFragments, // total number of reads we've looked at
        std::atomic<uint64_t>& numAssignedFragments, // total number of assigned reads
        std::atomic<uint64_t>& batchNum,
        bool initialRound,
        bool& burnedIn,
        double& logForgettingMass,
        std::mutex& ffMutex,
        FragmentLengthDistribution& fragLengthDist,
        uint32_t minMEMLength,
        uint32_t maxMEMOcc,
        double coverageThresh,
        bool greedyChain,
        std::mutex& iomutex,
        size_t numThreads) {

            std::vector<std::thread> threads;


            rl.checkValid();
            // If the read library is paired-end
            // ------ Paired-end --------
            if (rl.format().type == ReadType::PAIRED_END) {
                char* readFiles[] = { const_cast<char*>(rl.mates1().front().c_str()),
                    const_cast<char*>(rl.mates2().front().c_str()) };

                size_t maxReadGroup{1000}; // Number of files to read simultaneously
                size_t concurrentFile{2}; // Number of reads in each "job"
                paired_parser parser(4 * numThreads, maxReadGroup, concurrentFile,
                        readFiles, readFiles + 2);

                for(int i = 0; i < numThreads; ++i)  {
                    if (greedyChain) {
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<paired_parser, TranscriptHitList>(
                                    &parser,
                                    numObservedFragments,
                                    numAssignedFragments,
                                    idx,
                                    transcripts,
                                    batchNum,
                                    logForgettingMass,
                                    ffMutex,
                                    clusterForest,
                                    fragLengthDist,
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    iomutex,
                                    initialRound,
                                    burnedIn);
                        };
                        threads.emplace_back(threadFun);
                    } else {
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<paired_parser, FragmentList>(
                                    &parser,
                                    std::ref(numObservedFragments),
                                    std::ref(numAssignedFragments),
                                    idx,
                                    std::ref(transcripts),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffMutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex),
                                    initialRound,
                                    std::ref(burnedIn));
                        };
                        threads.emplace_back(threadFun);
                    }
                }

                for(int i = 0; i < numThreads; ++i)
                    threads[i].join();

            } // ------ Single-end --------
            else if (rl.format().type == ReadType::SINGLE_END) {

                char* readFiles[] = { const_cast<char*>(rl.unmated().front().c_str()) };
                size_t maxReadGroup{1000}; // Number of files to read simultaneously
                size_t concurrentFile{1}; // Number of reads in each "job"
                stream_manager streams( rl.unmated().begin(),
                        rl.unmated().end(), concurrentFile);

                single_parser parser(4 * numThreads, maxReadGroup, concurrentFile,
                        streams);

                for(int i = 0; i < numThreads; ++i)  {
                    if (greedyChain) {
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<single_parser, TranscriptHitList>( &parser,
                                    std::ref(numObservedFragments), std::ref(numAssignedFragments),
                                    idx,
                                    std::ref(transcripts),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffMutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex),
                                    initialRound,
                                    std::ref(burnedIn));
                        };
                        threads.emplace_back(threadFun);
                    } else {
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<single_parser, FragmentList>( &parser,
                                    std::ref(numObservedFragments), std::ref(numAssignedFragments),
                                    idx,
                                    std::ref(transcripts),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffMutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex),
                                    initialRound,
                                    std::ref(burnedIn));
                        };
                        threads.emplace_back(threadFun);
                    }
                }
                for(int i = 0; i < numThreads; ++i)
                    threads[i].join();
            } // ------ END Single-end --------
}



/**
  *  Quantify the targets given in the file `transcriptFile` using the
  *  reads in the given set of `readLibraries`, and write the results
  *  to the file `outputFile`.  The reads are assumed to be in the format
  *  specified by `libFmt`.
  *
  */
void quantifyLibrary(
        ReadExperiment& experiment,
        bool greedyChain,
        uint32_t minMEMLength,
        uint32_t maxMEMOcc,
        double coverageThresh,
        size_t numRequiredFragments,
        uint32_t numQuantThreads) {

    bool burnedIn{false};
    //ErrorModel errMod(1.00);
    auto& refs = experiment.transcripts();
    size_t numTranscripts = refs.size();
    size_t miniBatchSize{1000};
    std::atomic<uint64_t> numObservedFragments{0};

    size_t maxFragLen = 800;
    size_t meanFragLen = 200;
    size_t fragLenStd = 80;
    size_t fragLenKernelN = 4;
    double fragLenKernelP = 0.5;
    FragmentLengthDistribution fragLengthDist(1.0, maxFragLen,
                                              meanFragLen, fragLenStd,
                                              fragLenKernelN,
                                              fragLenKernelP, 1);
    double logForgettingMass{std::log(1.0)};
    double forgettingFactor{0.60};
    bool initialRound{true};

    std::mutex ffMutex;
    std::mutex ioMutex;


    while (numObservedFragments < numRequiredFragments) {
        if (!initialRound) {
            if (!experiment.reset()) {
                fmt::print(stderr,
                  "\n\n======== WARNING ========\n"
                  "One of the provided read files: [{}] "
                  "is not a regular file and therefore can't be read from "
                  "more than once.\n\n"
                  "We observed only {} mapping fragments when we wanted at least {}.\n\n"
                  "Please consider re-running Salmon with these reads "
                  "as a regular file!\n"
                  "==========================\n\n",
                  experiment.readFilesAsString(), numObservedFragments, numRequiredFragments);
                break;
            }
        }

        auto processReadLibraryCallback =  [&](
                ReadLibrary& rl, bwaidx_t* idx,
                std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
                std::atomic<uint64_t>& numAssignedFragments,
                std::atomic<uint64_t>& batchNum, size_t numQuantThreads,
                bool& burnedIn) -> void  {
                processReadLibrary(rl, idx, transcripts, clusterForest,
                                   numObservedFragments, numAssignedFragments, batchNum,
                                   initialRound, burnedIn, logForgettingMass, ffMutex, fragLengthDist,
                                   minMEMLength, maxMEMOcc, coverageThresh, greedyChain,
                                   ioMutex, numQuantThreads);
              };

        // Process all of the reads
        experiment.processReads(numQuantThreads, processReadLibraryCallback);

        initialRound = false;
        //numObservedFragments += experiment.numObservedFragments();
        fmt::print(stderr, "\n# observed = {} / # required = {}\033[F",
                   numObservedFragments, numRequiredFragments);
    }
    fmt::print(stderr, "\n\n");
}

int salmonQuantify(int argc, char *argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool noBiasCorrect{false};
    bool optChain{false};
    uint32_t maxThreads = std::thread::hardware_concurrency();
    uint32_t sampleRate{1};
    size_t requiredObservations;
    double coverageThresh;
    uint32_t minMEMLength, maxMEMOcc;
    vector<string> unmatedReadFiles;
    vector<string> mate1ReadFiles;
    vector<string> mate2ReadFiles;

    po::options_description generic("salmon quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "sailfish index.")
    ("libtype,l", po::value<std::string>(), "Format string describing the library type")
    ("unmated_reads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
     "List of files containing unmated reads of (e.g. single-end reads)")
    ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
        "File containing the #1 mates")
    ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
        "File containing the #2 mates")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use concurrently.")
    //("sample,s", po::value<uint32_t>(&sampleRate)->default_value(1), "Sample rate --- only consider every s-th k-mer in a read.")
    ("num_required_obs,n", po::value(&requiredObservations)->default_value(50000000),
                                        "The minimum number of observations (mapped reads) that must be observed before\n"
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the \n"
                                        "input file, then it will be read through multiple times.")
    ("minLen,k", po::value<uint32_t>(&minMEMLength)->default_value(15), "MEMs smaller than this size won't be considered")
    ("maxOcc,m", po::value<uint32_t>(&maxMEMOcc)->default_value(100), "MEMs occuring more than this many times won't be considered")
    ("coverage,c", po::value<double>(&coverageThresh)->default_value(0.75), "required coverage of read by union of MEMs to consider it a \"hit\"")
    ("output,o", po::value<std::string>(), "Output quantification file")
    ("optChain", po::bool_switch(&optChain)->default_value(false), "Chain MEMs optimally rather than greed")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
    ("gene_map,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided\n"
                                        "Sailfish will output both quant.sf and quant.genes.sf files, where the latter\n"
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping\n"
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format\n"
                                        "where each line contains the name of a transcript and the gene to which it belongs\n"
                                        "separated by a tab.  The extension of the file is used to determine how the file\n"
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF\n"
                                        "format; files with any other extension are assumed to be in the simple format");

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(generic).run();

        po::store(orderedOptions, vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Quant
==========
Perform streaming SMEM-based estimation of
transcript abundance from RNA-seq reads
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }

        po::notify(vm);

        std::stringstream commentStream;
        commentStream << "# salmon (smem-based) v" << salmon::version << "\n";
        commentStream << "# [ program ] => salmon \n";
        commentStream << "# [ command ] => quant \n";
        for (auto& opt : orderedOptions.options) {
            commentStream << "# [ " << opt.string_key << " ] => {";
            for (auto& val : opt.value) {
                commentStream << " " << val;
            }
            commentStream << " }\n";
        }
        std::string commentString = commentStream.str();
        fmt::print(stderr, "{}", commentString);

        // Verify the gene_map before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("gene_map")) {
            // Make sure the provided file exists
            geneMapPath = vm["gene_map"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                std::cerr << "Could not fine transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'gene_map\' option or provide a valid file\n";
                std::exit(1);
            }
        }


        bool greedyChain = !optChain;
        bfs::path outputDirectory(vm["output"].as<std::string>());
        bfs::create_directory(outputDirectory);
        if (!(bfs::exists(outputDirectory) and bfs::is_directory(outputDirectory))) {
            std::cerr << "Couldn't create output directory " << outputDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }

        bfs::path indexDirectory(vm["index"].as<string>());
        bfs::path logDirectory = outputDirectory.parent_path() / "logs";

#if HAVE_LOGGER
        bfs::create_directory(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "writing logs to " << logDirectory.string() << "\n";
        g2LogWorker logger(argv[0], logDirectory.string());
        g2::initializeLogging(&logger);
#endif

        LOG(INFO) << "parsing read library format";

        vector<ReadLibrary> readLibraries = sailfish::utils::extractReadLibraries(orderedOptions);
        ReadExperiment experiment(readLibraries, indexDirectory);
        uint32_t nbThreads = vm["threads"].as<uint32_t>();

        quantifyLibrary(experiment, greedyChain, minMEMLength, maxMEMOcc, coverageThresh,
                        requiredObservations, nbThreads);

        // ---- Get rid of things we no longer need --------
        // bwa_idx_destroy(idx);

        size_t tnum{0};

        std::cerr << "writing output \n";

        bfs::path estFilePath = outputDirectory / "quant.sf";
        salmon::utils::writeAbundances(experiment, estFilePath, commentString);

        // Not ready for read / alignment-based bias correction yet
        if (!noBiasCorrect) {
            fmt::print(stderr, "Post-hoc bias correction is not yet supported in salmon; disabling\n");
            noBiasCorrect = true;
        }

        if (!noBiasCorrect) {
            // Estimated read length
            double estimatedReadLength = 36;
            // Number of k-mers per read
            double kmersPerRead = 1.0;
            // Total number of mapped kmers
            uint64_t mappedKmers= experiment.numMappedReads();

            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasFeatPath = indexDirectory / "bias_feats.txt";
            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrection(biasFeatPath, estFilePath, estimatedReadLength, kmersPerRead, mappedKmers,
                    20, biasCorrectedFile, nbThreads);

        }

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("gene_map")) {
            try {
                sailfish::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory,
                                                            !noBiasCorrect);
            } catch (std::invalid_argument& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                                   "estimates. The gene-level file(s) may not exist",
                                   e.what());
            }
        }

    } catch (po::error &e) {
        std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
    }


    return 0;
}

/*
int salmonQuantifyOrig(int argc, char *argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool noBiasCorrect{false};
    bool optChain{false};
    uint32_t maxThreads = std::thread::hardware_concurrency();
    uint32_t sampleRate{1};
    double coverageThresh;
    uint32_t minMEMLength, maxMEMOcc;
    vector<string> unmatedReadFiles;
    vector<string> mate1ReadFiles;
    vector<string> mate2ReadFiles;

    po::options_description generic("salmon quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>(), "sailfish index.")
    ("libtype,l", po::value<std::string>(), "Format string describing the library type")
    ("unmated_reads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
     "List of files containing unmated reads of (e.g. single-end reads)")
    ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
        "File containing the #1 mates")
    ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
        "File containing the #2 mates")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use concurrently.")
    //("sample,s", po::value<uint32_t>(&sampleRate)->default_value(1), "Sample rate --- only consider every s-th k-mer in a read.")
    ("minLen,k", po::value<uint32_t>(&minMEMLength)->default_value(15), "MEMs smaller than this size won't be considered")
    ("maxOcc,n", po::value<uint32_t>(&maxMEMOcc)->default_value(100), "MEMs occuring more than this many times won't be considered")
    ("coverage,c", po::value<double>(&coverageThresh)->default_value(0.75), "required coverage of read by union of MEMs to consider it a \"hit\"")
    ("output,o", po::value<std::string>(), "Output quantification file")
    ("optChain", po::bool_switch(&optChain)->default_value(false), "Chain MEMs optimally rather than greed")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
    ("gene_map,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided\n"
                                        "Sailfish will output both quant.sf and quant.genes.sf files, where the latter\n"
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping\n"
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format\n"
                                        "where each line contains the name of a transcript and the gene to which it belongs\n"
                                        "separated by a tab.  The extension of the file is used to determine how the file\n"
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF\n"
                                        "format; files with any other extension are assumed to be in the simple format");

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(generic).run();

        po::store(orderedOptions, vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Quant
==========
Perform streaming SMEM-based estimation of
transcript abundance from RNA-seq reads
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }

        po::notify(vm);


        for (auto& opt : orderedOptions.options) {
            std::cerr << "[ " << opt.string_key << " ] => {";
            for (auto& val : opt.value) {
                std::cerr << " " << val;
            }
            std::cerr << " }\n";
        }

        // Verify the gene_map before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("gene_map")) {
            // Make sure the provided file exists
            geneMapPath = vm["gene_map"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                std::cerr << "Could not fine transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'gene_map\' option or provide a valid file\n";
                std::exit(1);
            }
        }


        bool greedyChain = !optChain;
        bfs::path outputDirectory(vm["output"].as<std::string>());
        bfs::create_directory(outputDirectory);
        if (!(bfs::exists(outputDirectory) and bfs::is_directory(outputDirectory))) {
            std::cerr << "Couldn't create output directory " << outputDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }

        bfs::path indexDirectory(vm["index"].as<string>());
        bfs::path logDirectory = outputDirectory.parent_path() / "logs";

#if HAVE_LOGGER
        bfs::create_directory(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "writing logs to " << logDirectory.string() << "\n";
        g2LogWorker logger(argv[0], logDirectory.string());
        g2::initializeLogging(&logger);
#endif

        LOG(INFO) << "parsing read library format";

        vector<ReadLibrary> readLibraries = sailfish::utils::extractReadLibraries(orderedOptions);

        bwaidx_t *idx{nullptr};

        { // mem-based
            bfs::path indexPath = indexDirectory / "bwaidx";
            if ((idx = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS)) == 0) {
                fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                return 1;
            }
        }

        size_t numRecords = idx->bns->n_seqs;
        std::vector<Transcript> transcripts_tmp;

        fmt::print(stderr, "Index contained {} targets\n", numRecords);
        //transcripts_.resize(numRecords);
        for (auto i : boost::irange(size_t(0), numRecords)) {
            uint32_t id = i;
            char* name = idx->bns->anns[i].name;
            uint32_t len = idx->bns->anns[i].len;
            // copy over the length, then we're done.
            transcripts_tmp.emplace_back(id, name, len);
        }

        std::sort(transcripts_tmp.begin(), transcripts_tmp.end(),
                [](const Transcript& t1, const Transcript& t2) -> bool {
                    return t1.id < t2.id;
                });

        std::vector<Transcript> transcripts_;
        for (auto& t : transcripts_tmp) {
            transcripts_.emplace_back(t.id, t.RefName.c_str(), t.RefLength);
        }
        transcripts_tmp.clear();
        // --- done ---

        std::vector<std::thread> threads;
        std::atomic<uint64_t> totalHits(0);
        std::atomic<uint64_t> rn{0};
        std::atomic<uint64_t> processedReads{0};
        std::atomic<uint64_t> batchNum{0};

        size_t numTranscripts{transcripts_.size()};
        ClusterForest clusterForest(numTranscripts, transcripts_);

        size_t maxFragLen = 800;
        size_t meanFragLen = 200;
        size_t fragLenStd = 80;
        size_t fragLenKernelN = 4;
        double fragLenKernelP = 0.5;
        FragmentLengthDistribution fragLengthDist(1.0, maxFragLen, meanFragLen, fragLenStd, fragLenKernelN, fragLenKernelP, 1);

        double logForgettingMass{sailfish::math::LOG_1};
        std::mutex ffmutex;
        std::mutex iomutex;

        uint32_t nbThreads = vm["threads"].as<uint32_t>();
        for (auto& rl : readLibraries) {

            rl.checkValid();
            // If the read library is paired-end
            // ------ Paired-end --------
            if (rl.format().type == ReadType::PAIRED_END) {
                char* readFiles[] = { const_cast<char*>(rl.mates1().front().c_str()),
                    const_cast<char*>(rl.mates2().front().c_str()) };

                size_t maxReadGroup{1000}; // Number of files to read simultaneously
                size_t concurrentFile{2}; // Number of reads in each "job"
                paired_parser parser(4 * nbThreads, maxReadGroup, concurrentFile,
                        readFiles, readFiles + 2);


                for(int i = 0; i < nbThreads; ++i)  {
                    if (greedyChain) {
                        threads.push_back(std::thread(processReadsMEM<paired_parser, TranscriptHitList>, &parser,
                                    std::ref(totalHits), std::ref(rn), std::ref(processedReads),
                                    idx,
                                    std::ref(transcripts_),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffmutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex)
                                    ));
                    } else {
                        threads.push_back(std::thread(processReadsMEM<paired_parser, FragmentList>, &parser,
                                    std::ref(totalHits), std::ref(rn), std::ref(processedReads),
                                    idx,
                                    std::ref(transcripts_),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffmutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex)
                                    ));
                    }
                }

                for(int i = 0; i < nbThreads; ++i)
                    threads[i].join();

            } // ------ Single-end --------
            else if (rl.format().type == ReadType::SINGLE_END) {

                char* readFiles[] = { const_cast<char*>(rl.unmated().front().c_str()) };
                size_t maxReadGroup{1000}; // Number of files to read simultaneously
                size_t concurrentFile{1}; // Number of reads in each "job"
                stream_manager streams( rl.unmated().begin(),
                        rl.unmated().end(), concurrentFile);

                single_parser parser(4 * nbThreads, maxReadGroup, concurrentFile,
                        streams);

                for(int i = 0; i < nbThreads; ++i)  {
                    if (greedyChain) {
                        threads.push_back(std::thread(processReadsMEM<single_parser, TranscriptHitList>, &parser,
                                    std::ref(totalHits), std::ref(rn), std::ref(processedReads),
                                    idx,
                                    std::ref(transcripts_),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffmutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex)
                                    ));
                    } else {
                        threads.push_back(std::thread(processReadsMEM<single_parser, FragmentList>, &parser,
                                    std::ref(totalHits), std::ref(rn), std::ref(processedReads),
                                    idx,
                                    std::ref(transcripts_),
                                    std::ref(batchNum),
                                    std::ref(logForgettingMass),
                                    std::ref(ffmutex),
                                    std::ref(clusterForest),
                                    std::ref(fragLengthDist),
                                    minMEMLength,
                                    maxMEMOcc,
                                    coverageThresh,
                                    std::ref(iomutex)
                                    ));
                    }

                }

                for(int i = 0; i < nbThreads; ++i)
                    threads[i].join();
            } // ------ END Single-end --------

            std::cerr << "\n\n";
            std::cerr << "processed " << rn << " total reads\n";
            std::cout << "Had a hit for " << totalHits  / static_cast<double>(rn) * 100.0 << "% of the reads\n";

        }

        // ---- Get rid of things we no longer need --------
        bwa_idx_destroy(idx);

        size_t tnum{0};

        std::cerr << "writing output \n";

        bfs::path estFilePath = outputDirectory / "quant.sf";

        {
            std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(estFilePath.c_str(), "w"), std::fclose);

            fmt::print(output.get(), "# Salmon version {}\n", salmon::version);
            fmt::print(output.get(), "# Name\tLength\tTPM\tFPKM\tNumReads\n");

            const double million = 1000000.0;
            const double logBillion = std::log(1000000000.0);
            const double logNumFragments = std::log(static_cast<double>(rn));
            auto clusters = clusterForest.getClusters();
            size_t clusterID = 0;
            for(auto cptr : clusters) {
                double logClusterMass = cptr->logMass();
                double logClusterCount = std::log(static_cast<double>(cptr->numHits()));

                if (logClusterMass == sailfish::math::LOG_0) {
                    std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
                }

                bool requiresProjection{false};

                auto& members = cptr->members();
                size_t clusterSize{0};
                for (auto transcriptID : members) {
                    Transcript& t = transcripts_[transcriptID];
                    t.uniqueCounts = t.uniqueCount();
                    t.totalCounts = t.totalCount();
                }

                for (auto transcriptID : members) {
                    Transcript& t = transcripts_[transcriptID];
                    double logTranscriptMass = t.mass();
                    double logClusterFraction = logTranscriptMass - logClusterMass;
                    t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
                    requiresProjection |= t.projectedCounts > static_cast<double>(t.totalCounts) or
                        t.projectedCounts < static_cast<double>(t.uniqueCounts);
                    ++clusterSize;
                }

                if (clusterSize > 1 and requiresProjection) {
                    cptr->projectToPolytope(transcripts_);
                }

                ++clusterID;
            }

            double tfracDenom{0.0};
            double numMappedReads = rn;
            for (auto& transcript : transcripts_) {
                tfracDenom += (transcript.projectedCounts / numMappedReads) / transcript.RefLength;
            }

            // Now posterior has the transcript fraction
            for (auto& transcript : transcripts_) {
                double logLength = std::log(transcript.RefLength);
                double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
                double count = transcript.projectedCounts;
                //double countTotal = transcripts_[transcriptID].totalCounts;
                //double countUnique = transcripts_[transcriptID].uniqueCounts;
                double fpkm = count > 0 ? fpkmFactor * count : 0.0;
                double npm = (transcript.projectedCounts / numMappedReads);
                double tfrac = (npm / transcript.RefLength) / tfracDenom;
                double tpm = tfrac * million;

                fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
                        transcript.RefName, transcript.RefLength,
                        tpm, fpkm, count);
            }

        }

        // Not ready for read / alignment-based bias correction yet
        if (!noBiasCorrect) {
            fmt::print(stderr, "Post-hoc bias correction is not yet supported in salmon; disabling\n");
            noBiasCorrect = true;
        }

        if (!noBiasCorrect) {
            // Estimated read length
            double estimatedReadLength = 36;
            // Number of k-mers per read
            double kmersPerRead = 1.0;
            // Total number of mapped kmers
            uint64_t mappedKmers= totalHits;

            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasFeatPath = indexDirectory / "bias_feats.txt";
            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrection(biasFeatPath, estFilePath, estimatedReadLength, kmersPerRead, mappedKmers,
                    20, biasCorrectedFile, nbThreads);

        }

        // If the user requested gene-level abundances, then compute those now
        if (vm.count("gene_map")) {
            try {
                sailfish::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory,
                                                            !noBiasCorrect);
            } catch (std::invalid_argument& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                                   "estimates. The gene-level file(s) may not exist",
                                   e.what());
            }
        }

    } catch (po::error &e) {
        std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
    }


    return 0;
}
*/
