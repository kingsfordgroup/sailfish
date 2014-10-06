/**
>HEADER
    Copyright (c) 2013, 2014 Rob Patro rob.patro@cs.stonybrook.edu

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

// C++ string formatting library
#include "format.h"

// C Includes for BWA
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
}

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

// Boost Includes
#include <boost/filesystem.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>

// TBB Includes
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

// Logger includes
#if HAVE_LOGGER
#include "g2logworker.h"
#include "g2log.h"
#endif

// Cereal includes
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

// Sailfish / Salmon includes
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


class SMEMAlignment {
    public:
        SMEMAlignment(TranscriptID transcriptIDIn, LibraryFormat format,
                  double scoreIn = 0.0, uint32_t fragLengthIn= 0,
                  double logProbIn = sailfish::math::LOG_0) :
            transcriptID_(transcriptIDIn), format_(format), score_(scoreIn),
            fragLength_(fragLengthIn), logProb(logProbIn) {}

        inline TranscriptID transcriptID() { return transcriptID_; }
        inline uint32_t fragLength() { return fragLength_; }
        inline LibraryFormat libFormat() { return format_; }
        inline double score() { return score_; }
        // inline double coverage() {  return static_cast<double>(kmerCount) / fragLength_; };

        uint32_t kmerCount;
        double logProb;

    private:
        TranscriptID transcriptID_;
        LibraryFormat format_;
        double score_;
        uint32_t fragLength_;
};


/*
class FragAlignmentGroup {
    public:
        FragAlignmentGroup() : read_(nullptr) {}

        void setRead(std::string* r) { read_ = r; }
        std::string* read() { return read_; }
        std::vector<Alignment>& alignments() { return alignments_; }
        size_t numAlignments() { return alignments_.size(); }
        size_t size() { return alignments_.size(); }
        std::vector<Alignment>::iterator begin() { return alignments_.begin(); }
        std::vector<Alignment>::iterator end() { return alignments_.end(); }
    private:
        std::string* read_;
        std::vector<Alignment> alignments_;
};
*/


inline double logAlignFormatProb(const LibraryFormat observed, const LibraryFormat expected) {
    if (observed.type != expected.type or
        observed.orientation != expected.orientation ) {
        //std::cerr << "Expected format " << expected << ", but observed " << observed << "\n";
        return sailfish::math::LOG_0;
    } else {
        if (expected.strandedness == ReadStrandedness::U) {
            return sailfish::math::LOG_ONEHALF;
        } else {
            if (expected.strandedness == observed.strandedness) {
                return sailfish::math::LOG_1;
            } else {
                std::cerr << "Expected format " << expected << ", but observed " << observed << "\n";
                return sailfish::math::LOG_0;
            }
        }
    }

    fmt::print(stderr, "WARNING: logAlignFormatProb --- should not get here");
    return sailfish::math::LOG_0;
}

void processMiniBatch(
        double logForgettingMass,
        const ReadLibrary& readLib,
        std::vector<AlignmentGroup<SMEMAlignment>>& batchHits,
        //std::vector<std::vector<Alignment>>& batchHits,
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
    using sailfish::math::LOG_ONEHALF;
    using sailfish::math::logAdd;
    using sailfish::math::logSub;

    constexpr uint64_t numBurninFrags = 5000000;

    size_t numTranscripts{transcripts.size()};
    size_t localNumAssignedFragments{0};
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

    bool updateCounts = initialRound;

    const auto expectedLibraryFormat = readLib.format();

    // Build reverse map from transcriptID => hit id
    using HitID = uint32_t;
    btree::btree_map<TranscriptID, std::vector<SMEMAlignment*>> hitsForTranscript;
    size_t hitID{0};
    for (auto& hv : batchHits) {
        for (auto& tid : hv.alignments()) {
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
            // AGCHANGE: auto firstTranscriptID = alnGroup.front().transcriptID();
            auto firstTranscriptID = alnGroup.alignments().front().transcriptID();
            std::unordered_set<size_t> observedTranscripts;
            for (auto& aln : alnGroup.alignments()) {
                auto transcriptID = aln.transcriptID();
                auto& transcript = transcripts[transcriptID];
                transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;
                double coverage = aln.score();
                double logFragProb = (coverage > 0) ? std::log(coverage) : LOG_0;

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
                    // The probability that the fragments align to the given strands in the
                    // given orientations.
                    double logAlignCompatProb = logAlignFormatProb(aln.libFormat(), expectedLibraryFormat);

                    //aln.logProb = std::log(aln.kmerCount) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    //aln.logProb = std::log(std::pow(aln.kmerCount,2.0)) + (transcriptLogCount - logRefLength);// + qualProb + errLike;
                    aln.logProb = (transcriptLogCount - logRefLength) + logFragProb + logAlignCompatProb;// + qualProb + errLike;

                    sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);

                    if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                        if (updateCounts) { transcripts[transcriptID].addTotalCount(1); }
                        observedTranscripts.insert(transcriptID);
                    }
                } else {
                    aln.logProb = LOG_0;
                }
            }

            // If this fragment has a zero probability,
            // go to the next one
            if (sumOfAlignProbs == LOG_0) {
                continue;
            } else { // otherwise, count it as assigned
                ++localNumAssignedFragments;
            }

            // normalize the hits
            for (auto& aln : alnGroup.alignments()) {
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
                clusterForest.mergeClusters<SMEMAlignment>(alnGroup.alignments().begin(), alnGroup.alignments().end());
                clusterForest.updateCluster(alnGroup.alignments().front().transcriptID(), 1, logForgettingMass, updateCounts);
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
                std::for_each(hits.begin(), hits.end(), [&](SMEMAlignment* aln) -> void {
                        if (!std::isfinite(aln->logProb)) { std::cerr << "hitMass = " << aln->logProb << "\n"; }
                        hitMass = logAdd(hitMass, aln->logProb);
                        });

                double updateMass = logForgettingMass + hitMass;
                individualTotal = logAdd(individualTotal, updateMass);
                totalMass = logAdd(totalMass, updateMass);
                transcript.addMass(updateMass);
            } // end for
        } // end timer
        numAssignedFragments += localNumAssignedFragments;
        if (numAssignedFragments >= numBurninFrags and !burnedIn) { burnedIn = true; }
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
        KmerVote(int32_t vp, uint32_t rp, uint32_t vl) : votePos(vp), readPos(rp), voteLen(vl) {}
        int32_t votePos{0};
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
        int32_t bestHitPos{0};
        uint32_t bestHitCount{0};
        double bestHitScore{0.0};

        std::vector<KmerVote> votes;
        std::vector<KmerVote> rcVotes;

        uint32_t targetID;

        bool isForward_{true};

        void addFragMatch(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            int32_t votePos = static_cast<int32_t>(tpos) - readPos;
            votes.emplace_back(votePos, readPos, voteLen);
        }

        void addFragMatchRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            int32_t votePos = static_cast<int32_t>(tpos) - (readPos + 1);
            rcVotes.emplace_back(votePos, readPos, voteLen);
        }

        uint32_t totalNumHits() { return std::max(votes.size(), rcVotes.size()); }

        bool computeBestLoc_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                             std::string& read, bool isRC,
                             int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {
            // Did we update the highest-scoring cluster? This will be set to
            // true iff we have a cluster of a higher score than the score
            // currently given in maxClusterCount.
            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            struct VoteInfo {
                uint32_t coverage = 0;
                int32_t rightmostBase = 0;
            };

            uint32_t readLen = read.length();

            boost::container::flat_map<uint32_t, VoteInfo> hitMap;
            int32_t currClust{static_cast<int32_t>(sVotes.front().votePos)};
            for (size_t j = 0; j < sVotes.size(); ++j) {

                int32_t votePos = sVotes[j].votePos;
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
                    maxClusterScore = maxClusterCount / static_cast<double>(readLen);
                    updatedMaxScore = true;
                }

            }
            return updatedMaxScore;
        }

        bool computeBestLoc2_(std::vector<KmerVote>& sVotes, uint32_t tlen,
                              int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {

            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            double weights[] = { 1.0, 0.983471453822, 0.935506985032,
                0.860707976425, 0.765928338365, 0.6592406302, 0.548811636094,
                0.441902209585, 0.344153786865, 0.259240260646,
                0.188875602838};

            uint32_t maxGap = 4;
            uint32_t leftmost = (sVotes.front().votePos > maxGap) ? (sVotes.front().votePos - maxGap) : 0;
            uint32_t rightmost = std::min(sVotes.back().votePos + maxGap, tlen);

            uint32_t span = (rightmost - leftmost);
            std::vector<double> probAln(span, 0.0);
            double kwidth = 1.0 / (2.0 * maxGap);

            size_t nvotes = sVotes.size();
            for (size_t j = 0; j < nvotes; ++j) {
                uint32_t votePos = sVotes[j].votePos;
                uint32_t voteLen = sVotes[j].voteLen;

                auto x = j + 1;
                while (x < nvotes and sVotes[x].votePos == votePos) {
                    voteLen += sVotes[x].voteLen;
                    j += 1;
                    x += 1;
                }


                uint32_t dist{0};
                size_t start = (votePos >= maxGap) ? (votePos - maxGap - leftmost) : (votePos - leftmost);
                size_t mid = votePos - leftmost;
                size_t end = std::min(votePos + maxGap - leftmost, rightmost - leftmost);
                for (size_t k = start; k < end; k += 1) {
                    dist = (mid > k) ? mid - k : k - mid;
                    probAln[k] += weights[dist] * voteLen;
                    if (probAln[k] > maxClusterScore) {
                        maxClusterScore = probAln[k];
                        maxClusterPos = k + leftmost;
                        updatedMaxScore = true;
                    }
                }
            }

            return updatedMaxScore;
        }


        inline uint32_t numSampledHits_(Transcript& transcript, std::string& readIn,
                                        int32_t votePos, bool isRC, uint32_t numTries) {

            // The read starts at this position in the transcript (may be negative!)
            int32_t readStart = votePos;
            // The (uncorrected) length of the read
            int32_t readLen = readIn.length();
            // Pointer to the sequence of the read
            const char* read = readIn.c_str();
            // If the read starts before the first base of the transcript,
            // trim off the initial overhang  and correct the other variables
            if (readStart < 0) {
                uint32_t correction = -readStart;
                read += correction;
                readLen -= correction;
                readStart = 0;
            }
            // If the read hangs off the end of the transcript,
            // shorten its effective length.
            if (readStart + readLen >= transcript.RefLength) {
                readLen = transcript.RefLength - (readStart + 1);
            }

            // If the read is too short, it's not useful
            if (readLen <= 6) { return 0; }
            // The step between sample centers (given the number of samples we're going to take)
            double step = (readLen - 1) / static_cast<double>(numTries-1);
            // The strand of the transcript from which we'll extract sequence
            auto dir = (isRC) ? sailfish::stringtools::strand::reverse :
                                sailfish::stringtools::strand::forward;

            /*
            if (readStart > 100000) {
                std::cerr << "SUSPICIOUS READ START!!!!\n";
                std::cerr << "INFO: votePos = " << votePos << ", readStart = " << readStart << ", readLen = " << readLen << ", tran len = " << transcript.RefLength << "\n";
            }
            */

            /*
            std::stringstream ss;
            ss << "Supposed hit " << (isRC ? "RC" : "") << "\n";
            ss << "info: votePos = " << votePos << ", readStart = " << readStart << ", readLen = " << readLen << ", tran len = " << transcript.RefLength << ", step = " << step << "\n";
            if (readStart + readLen > transcript.RefLength) {
                ss << "ERROR!!!\n";
                std::cerr << "[[" << ss.str() << "]]";
                std::exit(1);
            } else {
                std::cerr << ss.str();
            }
            ss << "T : ";
            try {
            for ( size_t j = 0; j < readLen; ++j) {
                if (isRC) {
                    ss << transcript.charBaseAt(readStart+readLen-j,dir);
                } else {
                    ss << transcript.charBaseAt(readStart+j);
                }
            }
            ss << "\n";
            ss << "R : " << read << "\n";
            std::cerr << ss.str() << "\n";
            } catch (std::exception& e) {
                std::cerr << "EXCEPTION !!!!!! " << e.what() << "\n";
            }
            */

            // The index of the current sample within the read
            int32_t readIndex = 0;

            // The number of loci in the subvotes and their
            // offset patternns
            size_t lpos = 3;
            int leftPattern[] = {-4, -2, 0};
            int rightPattern[] = {0, 2, 4};
            int centerPattern[] = {-4, 0, 4};

            // The number of subvote hits we've had
            uint32_t numHits = 0;
            // Take the samples
            for (size_t i  = 0; i < numTries; ++i) {
                // The sample will be centered around this point
                readIndex = static_cast<uint32_t>(std::round(readStart + i * step)) - readStart;

                // The number of successful sub-ovtes we have
                uint32_t subHit = 0;
                // Select the center sub-vote pattern, unless we're near the end of a read
                int* pattern = &centerPattern[0];
                if (readIndex + pattern[0] < 0) {
                    pattern = &rightPattern[0];
                } else if (readIndex + pattern[lpos-1] >= readLen) {
                    pattern = &leftPattern[0];
                }

                // collect the subvotes
                for (size_t j = 0; j < lpos; ++j) {
                    // the pattern offset
                    int offset = pattern[j];
                    // and sample position it implies within the read
                    int readPos = readIndex + offset;

                    if (readStart + readPos >= transcript.RefLength) {
                        std::cerr  << "offset = " << offset << ", readPos = " << readPos << ", readStart = " << readStart << ", readStart + readPos = " << readStart + readPos << ", tlen = " << transcript.RefLength << "\n";
                    }

                    subHit += (isRC) ?
                        (transcript.charBaseAt(readStart + readLen - readPos, dir) == sailfish::stringtools::charCanon[read[readPos]]) :
                        (transcript.charBaseAt(readStart + readPos               ) == sailfish::stringtools::charCanon[read[readPos]]);
                }
                // if the entire subvote was successful, this is a hit
                numHits += (subHit == lpos);
            }
            // return the number of hits we had
            return numHits;
        }



        bool computeBestLoc3_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                              std::string& read, bool isRC,
                              int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {

            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            struct LocHitCount {
                int32_t loc;
                uint32_t nhits;
            };

            uint32_t numSamp = 15;
            std::vector<LocHitCount> hitCounts;
            size_t nvotes = sVotes.size();
            int32_t prevPos = -std::numeric_limits<int32_t>::max();
            for (size_t j = 0; j < nvotes; ++j) {
                int32_t votePos = sVotes[j].votePos;
                if (prevPos == votePos) { continue; }
                auto numHits = numSampledHits_(transcript, read, votePos, isRC, numSamp);
                hitCounts.push_back({votePos, numHits});
                prevPos = votePos;
            }

            uint32_t maxGap = 8;
            uint32_t hitIdx = 0;
            uint32_t accumHits = 0;
            int32_t hitLoc = hitCounts[hitIdx].loc;
            while (hitIdx < hitCounts.size()) {
                uint32_t idx2 = hitIdx;
                while (idx2 < hitCounts.size() and std::abs(hitCounts[idx2].loc - hitLoc) <= maxGap) {
                    accumHits += hitCounts[idx2].nhits;
                    ++idx2;
                }
                double score = static_cast<double>(accumHits) / numSamp;
                if (score > maxClusterScore) {
                    maxClusterCount = accumHits;
                    maxClusterScore = score;
                    maxClusterPos = hitCounts[hitIdx].loc;
                    updatedMaxScore = true;
                }
                accumHits = 0;
                ++hitIdx;
            }

            return updatedMaxScore;
        }


        bool computeBestChain(Transcript& transcript, std::string& read) {
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

            int32_t maxClusterPos{0};
            uint32_t maxClusterCount{0};
            double maxClusterScore{0.0};

            // we don't need the return value from the first call
            static_cast<void>(computeBestLoc3_(votes, transcript, read, false, maxClusterPos, maxClusterCount, maxClusterScore));
            bool revIsBest = computeBestLoc3_(rcVotes, transcript, read, true, maxClusterPos, maxClusterCount, maxClusterScore);
            isForward_ = not revIsBest;

            bestHitPos = maxClusterPos;
            bestHitCount = maxClusterCount;
            bestHitScore = maxClusterScore;
            //std::cerr << "SCORE: " << (isForward_ ? ("(FWD)") : ("(RC)")) << ", " << maxClusterScore << "\n";
            return true;
        }

        bool isForward() { return isForward_; }

        /*
        inline uint32_t numRandHits_(Transcript& t, std::string& read,
                                     bool isRC,
                                     uint32_t votePos, uint32_t readPos,
                                     uint32_t voteLen, uint32_t readLen,
                                     uint32_t numTries) {

            uint32_t transcriptOffset = votePos;
            uint32_t start1 = 0;
            uint32_t end1 = readPos;
            uint32_t start2 = readPos + voteLen;
            uint32_t end2 = start1 + readLen;
            uint32_t len1 = end1 - start1;
            uint32_t len2 = end2 - start2;
            uint32_t numSamp1 = 0;
            if (len2 == 0) {
                numSamp1 = numTries;
                numSamp2 = 0;
            }

            numSamp1 = std::round(numTries * (static_cast<double>(len1) / (len1+len2)));
            uint32_t step1 = std::floor(static_cast<double>(len1) / numSamp1);

            numSamp2 = numTries - numSamp1;
            uint32_t step2 = std::floor(static_cast<double>(len2) / numSamp2);

            uint32_t readIndex = start1;
            uint32_t transcriptIndex = start1 + offset;
            for (size_t i = 0; i < numSamp1; ++i) {
                numMatches += (isRC) ?
                              (transcript.charBaseAt(transcriptIndex) == sailfish::stringtools::charRC[read[readLen - readIndex]]) :
                              (transcript.charBaseAt(transcriptIndex) == read[readIndex]);
                readIndex += step1;
                transcriptIndex += step1;
            }

            readIndex = start2;
            transcriptIndex = start2 + offset;
            for (size_t i = 0; i < numSamp2; ++i) {
                numMatches += (isRC) ? () :
                              (transcript.charBaseAt(transcriptIndex) == sailfish::stringtools::charRC[read[readLen - readIndex]]) :
                              (transcript.charBaseAt(transcriptIndex) == read[readIndex]);
                readIndex += step2;
                transcriptIndex += step2;
            }

            return numMatches;
        }

        bool computeBestPosition_(std::vector<KmerVote>& sVotes,
                                  Transcript& target,
                                  std::string& read,
                                  bool isRC,
                                  uint32_t& maxClusterPos,
                                  uint32_t& maxClusterCount) {

            bool updatedMaxScore{false};
            if (sVotes.size() == 0) { return updatedMaxScore; }

            boost::container::flat_map<uint32_t, uint32_t> hitMap;
            uint32_t readLen = read.length();
            for (size_t j = 0; j < sVotes.size(); ++j) {
                uint32_t votePos = sVotes[j].votePos;
                uint32_t readPos = sVotes[j].readPos;
                uint32_t voteLen = sVotes[j].voteLen;
                uint32_t currCount = voteLen;

                if (currCount < readLen) {
                    currCount += numRandHits_(target, read, votePos, readPos, voteLen, readLen, 15);
                }
                hitMap[votePos] += currCount;
                auto newNumHits = hitMap[votePos];

                if (newNumHits > maxClusterCount) {
                    maxClusterCount = newNumHits;
                    maxClusterPos = votePos;
                    updatedMaxScore = true;
                }
            }
            return updatedMaxScore;
        }

        bool computeBestPosition(Transcript& target, std::string& read) {
            using sailfish::math::LOG_0;

            uint32_t maxHitPos{0};
            uint32_t maxHitCount{0};
            double maxHitScore{LOG_0};

            // we don't need the return value from the first call
            static_cast<void>(computeBestPosition_(votes, target, read, maxHitPos, maxHitCount));
            bool revIsBest = computeBestPosition_(rcVotes, target, read, maxHitPos, maxHitCount);
            isForward_ = not revIsBest;

            bestHitPos = maxClusterPos;
            bestHitCount = maxClusterCount;
            bestHitScore = -std::log(std::exp(0.25, bestHitCount));

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
                        AlignmentGroup<SMEMAlignment>& hitList,
                        //std::vector<Alignment>& hitList,
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

    //---------- End 1 ----------------------//
    {
        std::string readStr   = frag.first.seq;
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
                        // If the hit is to the reverse strand
                        if (isRev) {
                            pos -= len - 1;
                        }
                        bns_cnt_ambi(idx->bns, pos, len, &refID);
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset);

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

    //---------- End 2 ----------------------//
    {
        std::string readStr   = frag.second.seq;
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
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset);

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
    auto& alnList = hitList.alignments();
    alnList.clear();

    for (auto& tHitList : leftHits) {
        // Coverage score
        Transcript& t = transcripts[tHitList.first];
        tHitList.second.computeBestChain(t, frag.first.seq);
        ++leftHitCount;
    }

    double cutoffLeft{ coverageThresh };//* leftReadLength};
    double cutoffRight{ coverageThresh };//* rightReadLength};


    for (auto& tHitList : rightHits) {
        auto it = leftHits.find(tHitList.first);
        // Coverage score
        if (it != leftHits.end() and it->second.bestHitScore >= cutoffLeft) {
            Transcript& t = transcripts[tHitList.first];
            tHitList.second.computeBestChain(t, frag.second.seq);
            if (tHitList.second.bestHitScore < cutoffRight) { continue; }

            auto end1Start = it->second.bestHitPos;
            auto end2Start = tHitList.second.bestHitPos;

            double score = (it->second.bestHitScore + tHitList.second.bestHitScore) * 0.5;
            uint32_t fragLength = std::abs(static_cast<int32_t>(end1Start) -
                                           static_cast<int32_t>(end2Start + rightReadLength));

            bool end1IsForward = it->second.isForward();
            bool end2IsForward = tHitList.second.isForward();

            auto fmt = salmon::utils::hitType(it->second.bestHitPos, end1IsForward, tHitList.second.bestHitPos, end2IsForward);

            alnList.emplace_back(tHitList.first, fmt, score, fragLength);
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
                        AlignmentGroup<SMEMAlignment>& hitList,
                        //std::vector<Alignment>& hitList,
                        uint64_t& hitListCount,
                        std::vector<Transcript>& transcripts) {

    uint64_t leftHitCount{0};

    //std::unordered_map<uint64_t, TranscriptHitList> hits;
    std::unordered_map<uint64_t, CoverageCalculator> hits;

    uint32_t readLength{0};

    //---------- get hits ----------------------//
    {
        std::string readStr   = frag.seq;
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
                        long hitLoc = static_cast<long>(pos - idx->bns->anns[refID].offset);
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
    auto& alnList = hitList.alignments();
    alnList.clear();

    double cutoff{ coverageThresh };//* readLength};
    for (auto& tHitList : hits) {
        // Coverage score
        Transcript& t = transcripts[tHitList.first];
        tHitList.second.computeBestChain(t, frag.seq);
        if (tHitList.second.bestHitScore >= cutoff) {
            double score = tHitList.second.bestHitScore;
            bool isForward = tHitList.second.isForward();

            auto fmt = salmon::utils::hitType(tHitList.second.bestHitPos, isForward);

            alnList.emplace_back(tHitList.first, fmt, score);
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
               const ReadLibrary& rl,
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

  std::vector<AlignmentGroup<SMEMAlignment>> hitLists;
  //std::vector<std::vector<Alignment>> hitLists;
  hitLists.resize(5000);

  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  // Super-MEM iterator
  smem_i *itr = smem_itr_init(idx->bwt);
  const bwtintv_v *a = nullptr;
  int minLen{static_cast<int>(minMEMLength)};
  int minIWidth{static_cast<int>(maxMEMOcc)};
  int splitWidth{0};

  auto expectedLibType = rl.format();

  size_t locRead{0};
  while(true) {
    typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    hitLists.resize(j->nb_filled);
    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch

        //hitLists[i].setRead(&j->data[i]);

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
        if (hitList.size() > 100) { hitList.alignments().clear(); }

    } // end for i < j->nb_filled

    auto oldBatchNum = batchNum++;
    if (oldBatchNum > 1) {
        ffMutex.lock();
        logForgettingMass += forgettingFactor * std::log(static_cast<double>(oldBatchNum-1)) -
        std::log(std::pow(static_cast<double>(oldBatchNum), forgettingFactor) - 1);
        ffMutex.unlock();
    }

    processMiniBatch(logForgettingMass, rl, hitLists, transcripts, clusterForest,
                     fragLengthDist, numAssignedFragments, eng, initialRound, burnedIn);

    // At this point, the parser can re-claim the strings
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
                                    rl,
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
                        /*
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<paired_parser, FragmentList>(
                                    &parser,
                                    rl,
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
                        */
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
                                    rl,
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
                        /*
                        auto threadFun = [&]() -> void {
                                    processReadsMEM<single_parser, FragmentList>( &parser,
                                    rl,
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
                        */
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

    size_t numPrevObservedFragments = 0;

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
            numPrevObservedFragments = numObservedFragments;
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
        fmt::print(stderr, "\n# observed = {} / # required = {}\n",
                   numObservedFragments, numRequiredFragments);
        fmt::print(stderr, "# assigned = {} / # observed (this round) = {}\033[F\033[F",
                   experiment.numAssignedFragments(),
                   numObservedFragments - numPrevObservedFragments);
    }
    fmt::print(stderr, "\n\n\n\n");
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
    ("index,i", po::value<string>(), "Salmon index")
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
                                        "The minimum number of observations (mapped reads) that must be observed before "
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the "
                                        "input file, then it will be read through multiple times.")
    ("minLen,k", po::value<uint32_t>(&minMEMLength)->default_value(15), "SMEMs smaller than this size won't be considered")
    ("maxOcc,m", po::value<uint32_t>(&maxMEMOcc)->default_value(100), "SMEMs occuring more than this many times won't be considered")
    ("coverage,c", po::value<double>(&coverageThresh)->default_value(0.75), "required coverage of read by union of SMEMs to consider it a \"hit\"")
    ("output,o", po::value<std::string>(), "Output quantification file")
    ("optChain", po::bool_switch(&optChain)->default_value(false), "Chain MEMs optimally rather than greed")
    ("no_bias_correct", po::value(&noBiasCorrect)->zero_tokens(), "turn off bias correction")
    ("gene_map,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
                                        "Salmon will output both quant.sf and quant.genes.sf files, where the latter "
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format "
                                        "where each line contains the name of a transcript and the gene to which it belongs "
                                        "separated by a tab.  The extension of the file is used to determine how the file "
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
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

