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
#include "btree_map.h"
#include "btree_set.h"


/** BWA Includes */
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cctype>

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "kseq.h"
#include "utils.h"
}

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

#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "ClusterForest.hpp"
#include "PerfectHashIndex.hpp"
#include "LookUpTableUtils.hpp"
#include "SailfishMath.hpp"
#include "Transcript.hpp"
#include "LibraryFormat.hpp"
#include "SailfishUtils.hpp"
#include "ReadLibrary.hpp"

#include "PairSequenceParser.hpp"

//#include "google/dense_hash_map"

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
        Alignment(TranscriptID transcriptIDIn, uint32_t kCountIn = 1, double logProbIn = sailfish::math::LOG_0) :
            transcriptID_(transcriptIDIn), kmerCount(kCountIn), logProb(logProbIn) {}

        inline TranscriptID transcriptID() { return transcriptID_; }
        uint32_t kmerCount;
        double logProb;

    private:
        TranscriptID transcriptID_;
};

void processMiniBatch(
        double logForgettingMass,
        std::vector<std::vector<Alignment>>& batchHits,
        std::vector<Transcript>& transcripts,
        ClusterForest& clusterForest
        ) {

    using sailfish::math::LOG_0;
    using sailfish::math::logAdd;
    using sailfish::math::logSub;

    size_t numTranscripts{transcripts.size()};

    bool burnedIn{true};
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

                    //aln.logProb = std::log(aln.kmerCount) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    //aln.logProb = std::log(std::pow(aln.kmerCount,2.0)) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    aln.logProb = (transcriptLogCount - logRefLength);// + qualProb + errLike;


                    sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);
                    if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                        transcripts[transcriptID].addTotalCount(1);
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
                /*
                double r = uni(eng);
                if (!burnedIn and r < std::exp(aln.logProb)) {
                    errMod.update(aln, transcript, aln.logProb, logForgettingMass);
                    if (aln.fragType() == ReadType::PAIRED_END) {
                        double fragLength = aln.fragLen();//std::abs(aln.read1->core.pos - aln.read2->core.pos) + aln.read2->core.l_qseq;
                        fragLengthDist.addVal(fragLength, logForgettingMass);
                    }
                }
                */

            } // end normalize

            // update the single target transcript
            if (transcriptUnique) {
                transcripts[firstTranscriptID].addUniqueCount(1);
                clusterForest.updateCluster(firstTranscriptID, 1, logForgettingMass);
            } else { // or the appropriate clusters
                clusterForest.mergeClusters<Alignment>(alnGroup.begin(), alnGroup.end());
                clusterForest.updateCluster(alnGroup.front().transcriptID(), 1, logForgettingMass);
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
        //if (processedReads >= 5000000 and !burnedIn) { burnedIn = true; }
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


class TranscriptHitList {
    public:
        uint32_t bestHitPos{0};
        uint32_t bestHitScore{0};
        std::vector<KmerVote> votes;
        std::vector<KmerVote> rcVotes;

        void addVote(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            uint32_t votePos = (readPos > tpos) ? 0 : tpos - readPos;
            votes.emplace_back(votePos, readPos, voteLen);
        }

        void addVoteRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
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
            int32_t currClust{sVotes.front().votePos};
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

        /*

        void computeBestChain_(std::vector<KmerVote>& sVotes, uint32_t& maxClusterPos, uint32_t& maxClusterCount) {
            if (sVotes.size() == 0) { return; }
            std::vector<double> scores(sVotes.size(), -std::numeric_limits<double>::max());

            for (size_t i = 0; i < sVotes.size(); ++i) {
                auto& fragI = sVotes[i];
                for (size_t j = 0; j < i; ++j) {
                    auto& fragJ = sVotes[j];
                    if (precedes(fragJ, fragI)) {
                        auto score = scores[j] - gapCost(fragJ, fragI);
                        if (score > scores[i]) { scores[i] = score; }
                    }
                }
                scores[i] += weight(fragI);
                if (scores[i] > maxClusterCount) {
                    maxClusterCount = scores[i];
                    maxClusterPos = i;
                }
            }
        }

        bool computeBestChain() {
            return true;
        }

        */

        bool computeBestHit() {
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
};


void getHitsForFragment(std::pair<header_sequence_qual, header_sequence_qual>& frag,
                        bwaidx_t *idx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        int minLen,
                        int minIWidth,
                        int splitWidth,
                        double coverageThresh,
                        std::vector<Alignment>& hitList,
                        uint64_t& hitListCount) {

    uint64_t leftHitCount{0};

    std::unordered_map<uint64_t, TranscriptHitList> leftHits;
    std::unordered_map<uint64_t, TranscriptHitList> rightHits;

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
                long numHits = static_cast<long>(p->x[2]);

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
                            leftHits[refID].addVoteRC(hitLoc, qstart, len);
                        } else {
                            leftHits[refID].addVote(hitLoc, qstart, len);
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
                long numHits = static_cast<long>(p->x[2]);

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
                            rightHits[refID].addVoteRC(hitLoc, qstart, len);
                        } else {
                            rightHits[refID].addVote(hitLoc, qstart, len);
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
        tHitList.second.computeBestHit();
        ++leftHitCount;
    }

    double cutoffLeft{ coverageThresh * leftReadLength};
    double cutoffRight{ coverageThresh * rightReadLength};

    for (auto& tHitList : rightHits) {
        auto it = leftHits.find(tHitList.first);
        // Coverage score
        if (it != leftHits.end() and it->second.bestHitScore >= cutoffLeft) {
            tHitList.second.computeBestHit();
            if (tHitList.second.bestHitScore < cutoffRight) { continue; }
            uint32_t score = it->second.bestHitScore + tHitList.second.bestHitScore;
            hitList.emplace_back(tHitList.first, score);
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
void getHitsForFragment(jellyfish::header_sequence_qual& frag,
                        bwaidx_t *idx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        int minLen,
                        int minIWidth,
                        int splitWidth,
                        double coverageThresh,
                        std::vector<Alignment>& hitList,
                        uint64_t& hitListCount) {

    uint64_t leftHitCount{0};

    std::unordered_map<uint64_t, TranscriptHitList> hits;

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
                long numHits = static_cast<long>(p->x[2]);

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
                            hits[refID].addVoteRC(hitLoc, qstart, len);
                        } else {
                            hits[refID].addVote(hitLoc, qstart, len);
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
        tHitList.second.computeBestHit();
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
template <typename ParserT>
void processReadsMEM(ParserT* parser,
               std::atomic<uint64_t>& totalHits,
               std::atomic<uint64_t>& rn,
               bwaidx_t *idx,
               std::vector<Transcript>& transcripts,
               std::atomic<uint64_t>& batchNum,
               double& logForgettingMass,
               std::mutex& ffMutex,
               ClusterForest& clusterForest,
               uint32_t minMEMLength,
               uint32_t maxMEMOcc,
               double coverageThresh
               ) {
  uint64_t count_fwd = 0, count_bwd = 0;

  double forgettingFactor{0.65};

  std::vector<std::vector<Alignment>> hitLists;
  hitLists.resize(5000);

  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  // Super-MEM iterator
  smem_i *itr = smem_itr_init(idx->bwt);
  const bwtintv_v *a;
  int minLen{minMEMLength};
  int minIWidth{maxMEMOcc};
  int splitWidth{0};

  size_t locRead{0};
  while(true) {
    typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;          // If got nothing, quit

    hitLists.resize(j->nb_filled);
    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got

    /*
        std::unordered_map<uint64_t, TranscriptHitList> leftFwdHits;
        std::unordered_map<uint64_t, TranscriptHitList> leftBwdHits;

        std::unordered_map<uint64_t, TranscriptHitList> rightFwdHits;
        std::unordered_map<uint64_t, TranscriptHitList> rightBwdHits;

        uint32_t totFwdLeft = 0;
        uint32_t totBwdLeft = 0;
        uint32_t totFwdRight = 0;
        uint32_t totBwdRight = 0;
        uint32_t readLength = 0;

        uint32_t leftReadLength{0};
        uint32_t rightReadLength{0};

        uint64_t prevEQ{std::numeric_limits<uint64_t>::max()};

        //---------- Left End ----------------------//
        {
            std::string& readStr   = j->data[i].first.seq;
            uint32_t readLen      = j->data[i].first.seq.size();

            leftReadLength = readLen;

            for (int p = 0; p < readLen; ++p) {
                readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
            }

            char* readPtr = const_cast<char*>(readStr.c_str());

            smem_set_query(itr, readLen, reinterpret_cast<uint8_t*>(readPtr));

            // while there are more matches on the query
            while ((a = smem_next(itr, minLen<<1, splitWidth)) != 0) {
                for (size_t mi = 0; mi < a->n; ++mi) {
                    bwtintv_t *p = &a->a[mi];
                    if (static_cast<uint32_t>(p->info) - (p->info>>32) < minLen) continue;
                    uint32_t qstart = static_cast<uint32_t>(p->info>>32);
                    uint32_t qend = static_cast<uint32_t>(p->info);
                    long numHits = static_cast<long>(p->x[2]);

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
                                leftFwdHits[refID].addVoteRC(hitLoc, qstart, len);
                            } else {
                                leftFwdHits[refID].addVote(hitLoc, qstart, len);
                            }
                        } // for k
                    } // if <= minIWidth

                } // for mi
            } // for all query matches
        }

        hitLists[i].clear();
        auto& hitList = hitLists[i];
        prevEQ = std::numeric_limits<uint64_t>::max();
        //---------- Right End ----------------------//
        {
            std::string& readStr   = j->data[i].second.seq;
            uint32_t readLen      = j->data[i].second.seq.size();

            rightReadLength = readLen;

            for (int p = 0; p < readLen; ++p) {
                readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
            }

            char* readPtr = const_cast<char*>(readStr.c_str());
            smem_set_query(itr, readLen, reinterpret_cast<uint8_t*>(readPtr));

            // while there are more matches on the query
            while ((a = smem_next(itr, minLen<<1, splitWidth)) != 0) {
                for (size_t mi = 0; mi < a->n; ++mi) {
                    bwtintv_t *p = &a->a[mi];
                    if (static_cast<uint32_t>(p->info) - (p->info>>32) < minLen) continue;
                    uint32_t qstart = static_cast<uint32_t>(p->info>>32);
                    uint32_t qend = static_cast<uint32_t>(p->info);
                    long numHits = static_cast<long>(p->x[2]);

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
                                rightFwdHits[refID].addVoteRC(hitLoc, qstart, len);
                            } else {
                                rightFwdHits[refID].addVote(hitLoc, qstart, len);
                            }
                        } // for k
                    } // if <= minIWidth

                } // for mi
            } // for all query matches

        } // end right

        size_t readHits{0};
        std::unordered_map<TranscriptID, uint32_t> leftHitCounts;
        std::unordered_map<TranscriptID, uint32_t> rightHitCounts;

        auto& leftHits = leftFwdHits;
        auto& rightHits = rightFwdHits;

        for (auto& tHitList : leftHits) {
            // Coverage score
            tHitList.second.computeBestHit();
            ++leftHitCount;
        }

        double cutoffLeft{ coverageThresh * leftReadLength};
        double cutoffRight{ coverageThresh * rightReadLength};

        for (auto& tHitList : rightHits) {
            auto it = leftHits.find(tHitList.first);
            // Coverage score
            if (it != leftHits.end() and it->second.bestHitScore >= cutoffLeft) {
                tHitList.second.computeBestHit();
                if (tHitList.second.bestHitScore < cutoffRight) { continue; }
                uint32_t score = it->second.bestHitScore + tHitList.second.bestHitScore;
                hitList.emplace_back(tHitList.first, score);
                readHits += score;
                ++hitListCount;
            }
        }

        totalHits += hitList.size() > 0;
    */


        getHitsForFragment(j->data[i], idx, itr, a,
                           minLen, minIWidth, splitWidth, coverageThresh,
                           hitLists[i], hitListCount);
        auto& hitList = hitLists[i];
        totalHits += hitList.size() > 0;


        locRead++;
        ++rn;
        if (rn % 50000 == 0) {
            std::cerr << "\r\rprocessed read "  << rn;
/*            std::cerr << "\n leftHits.size() " << leftHits.size()
                << ", rightHits.size() " << rightHits.size() << ", hit list of size = " << hitList.size() << "\n";
                std::cerr << "average leftHits = " << leftHitCount / static_cast<float>(locRead)
               << ", average hitList = " << hitListCount / static_cast<float>(locRead) << "\n";
               */
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

    processMiniBatch(logForgettingMass, hitLists, transcripts, clusterForest);
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



int salmonQuantify(int argc, char *argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool noBiasCorrect{false};
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
Perform streaming k-mer-group-based estimation of
transcript abundance from RNA-seq reads
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }

        po::notify(vm);


        for (auto& opt : orderedOptions.options) {
            std::cerr << "[ " << opt.string_key << "] => {";
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

        vector<ReadLibrary> readLibraries;
        for (auto& opt : orderedOptions.options) {
            if (opt.string_key == "libtype") {
                LibraryFormat libFmt = sailfish::utils::parseLibraryFormatString(opt.value[0]);
                if (libFmt.check()) {
                    std::cerr << libFmt << "\n";
                } else {
                    std::stringstream ss;
                    ss << libFmt << " is invalid!";
                    throw std::invalid_argument(ss.str());
                }
                readLibraries.emplace_back(libFmt);
            } else if (opt.string_key == "mates1") {
                readLibraries.back().addMates1(opt.value);
            } else if (opt.string_key == "mates2") {
                readLibraries.back().addMates2(opt.value);
            } else if (opt.string_key == "unmated_reads") {
                readLibraries.back().addUnmated(opt.value);
            }
        }

        bwaidx_t *idx;

        { // mem-based
            bfs::path indexPath = indexDirectory / "bwaidx";
            if ((idx = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
        }

        size_t numRecords = idx->bns->n_seqs;
        std::vector<Transcript> transcripts_tmp;

        std::cerr << "Transcript LUT contained " << numRecords << " records\n";
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
        std::atomic<uint64_t> batchNum{0};

        size_t numTranscripts{transcripts_.size()};
        ClusterForest clusterForest(numTranscripts, transcripts_);

        double logForgettingMass{sailfish::math::LOG_1};
        std::mutex ffmutex;


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
                    threads.push_back(std::thread(processReadsMEM<paired_parser>, &parser,
                                std::ref(totalHits), std::ref(rn),
                                idx,
                                std::ref(transcripts_),
                                std::ref(batchNum),
                                std::ref(logForgettingMass),
                                std::ref(ffmutex),
                                std::ref(clusterForest),
                                minMEMLength,
                                maxMEMOcc,
                                coverageThresh
                                ));
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
                    threads.push_back(std::thread(processReadsMEM<single_parser>, &parser,
                                std::ref(totalHits), std::ref(rn),
                                idx,
                                std::ref(transcripts_),
                                std::ref(batchNum),
                                std::ref(logForgettingMass),
                                std::ref(ffmutex),
                                std::ref(clusterForest),
                                minMEMLength,
                                maxMEMOcc,
                                coverageThresh
                                ));
                }

                for(int i = 0; i < nbThreads; ++i)
                    threads[i].join();
            } // ------ END Single-end --------

            std::cerr << "\n\n";
            std::cerr << "processed " << rn << " total reads\n";
            std::cout << "Had a hit for " << totalHits  / static_cast<double>(rn) * 100.0 << "% of the reads\n";

        }

        size_t tnum{0};

        std::cerr << "writing output \n";

        bfs::path estFilePath = outputDirectory / "quant.sf";
        std::ofstream output(estFilePath.string());
        output << "# Salmon v0.1.0\n";
        output << "# Name\tLength\tFPKM\tNumReads\n";

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

            // Now posterior has the transcript fraction
            size_t tidx = 0;
            for (auto transcriptID : members) {
                auto& transcript = transcripts_[transcriptID];
                double logLength = std::log(transcript.RefLength);
                double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
                double count = transcripts_[transcriptID].projectedCounts;
                double countTotal = transcripts_[transcriptID].totalCounts;
                double countUnique = transcripts_[transcriptID].uniqueCounts;
                double fpkm = count > 0 ? fpkmFactor * count : 0.0;
                output << transcript.RefName << '\t' << transcript.RefLength
                        << '\t' << fpkm << '\t' << fpkm
                        << '\t' << fpkm << '\t' << count <<  '\t' << count << '\n';
                ++tidx;
            }

            ++clusterID;
        }
        output.close();


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


        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("gene_map")) {
           std::cerr << "Computing gene-level abundance estimates\n";
           bfs::path gtfExtension(".gtf");
           auto extension = geneMapPath.extension();

           TranscriptGeneMap tranGeneMap;
           // parse the map as a GTF file
           if (extension == gtfExtension) {
               // Using the custom GTF Parser
               //auto features = GTFParser::readGTFFile<TranscriptGeneID>(geneMapPath.string());
               //tranGeneMap = sailfish::utils::transcriptToGeneMapFromFeatures(features);

               // Using libgff
                tranGeneMap = sailfish::utils::transcriptGeneMapFromGTF(geneMapPath.string(), "gene_id");
           } else { // parse the map as a simple format files
               std::ifstream tgfile(geneMapPath.string());
                tranGeneMap = sailfish::utils::readTranscriptToGeneMap(tgfile);
                tgfile.close();
           }

           std::cerr << "There were " << tranGeneMap.numTranscripts() << " transcripts mapping to "
                     << tranGeneMap.numGenes() << " genes\n";

           sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, estFilePath);
           /** Create a gene-level summary of the bias-corrected estimates as well if these exist **/
           if (!noBiasCorrect) {
                bfs::path biasCorrectEstFilePath(estFilePath.parent_path());
                biasCorrectEstFilePath /= "quant_bias_corrected.sf";
                sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, biasCorrectEstFilePath);
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


