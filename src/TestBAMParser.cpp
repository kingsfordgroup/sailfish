//#include "MultithreadedBAMParser.hpp"
//#include <boost/thread/thread.hpp>

extern "C" {
#include "htslib/sam.h"
//#include "sam.h"
}

#include "format.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/lockfree/queue.hpp>

#include <omp.h>
#include <forward_list>
#include <tbb/atomic.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <random>
#include <memory>
#include <exception>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <condition_variable>

#include <tbb/concurrent_queue.h>

#include <boost/container/flat_map.hpp>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/program_options.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/normal.hpp>

#include "AlignmentGroup.hpp"
#include "MiniBatchInfo.hpp"
#include "BAMQueue.hpp"
#include "SailfishMath.hpp"
#include "FASTAParser.hpp"
#include "LibraryFormat.hpp"
#include "Transcript.hpp"
#include "ReadPair.hpp"
#include "ErrorModel.hpp"
#include "FragmentLengthDistribution.hpp"

namespace bfs = boost::filesystem;
using sailfish::math::LOG_0;
using sailfish::math::LOG_1;
using sailfish::math::logAdd;
using sailfish::math::logSub;

template <typename FragT>
using AlignmentBatch = std::vector<FragT>;

template <typename FragT>
using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

using PriorAbundanceVector = std::vector<double>;
using PosteriorAbundanceVector = std::vector<double>;

struct RefSeq {
    RefSeq(char* name, uint32_t len) : RefName(name), RefLength(len) {}
    std::string RefName;
    uint32_t RefLength;
};

/**
 * Multiple each element in the vector `vec` by the factor `scale`.
 */
template <typename T>
void scaleBy(std::vector<T>& vec, T scale) {
    std::for_each(vec.begin(), vec.end(), [scale](T& ele)->void { ele *= scale; });
}

class TranscriptCluster {
    friend class ClusterForest;
public:
    TranscriptCluster() : members_(std::list<size_t>()), count_({0}), logMass_(LOG_0), active_(true) {}
    TranscriptCluster(size_t initialMember) : members_(std::list<size_t>(1,initialMember)), count_({0}),
                                              logMass_(LOG_0), active_(true) {
    }

    void incrementCount(size_t num) { count_ += num; }
    void addMass(double logNewMass) { logMass_ = logAdd(logMass_, logNewMass); }
    void merge(TranscriptCluster& other) {
        members_.splice(members_.begin(), other.members_);
        logMass_ = logAdd(logMass_, other.logMass_);
        count_ += other.count_;
    }

    std::list<size_t>& members() { return members_; }
    size_t numHits() { return count_.load(); }
    bool isActive() { return active_; }
    void deactivate() { active_ = false; }
    double logMass() { return logMass_; }

    // Adopted from https://github.com/adarob/eXpress/blob/master/src/targets.cpp
    void projectToPolytope(std::vector<Transcript>& allTranscripts) {
        using sailfish::math::approxEqual;
        constexpr size_t maxIter = 1000;
        size_t iter{0};

        // The transcript belonging to this cluster
        double clusterCounts{static_cast<double>(count_)};
        std::vector<Transcript*> transcripts;
        for (auto tid : members_) {
            transcripts.push_back(&allTranscripts[tid]);
        }
        // The cluster size
        size_t clusterSize = transcripts.size();

        boost::dynamic_bitset<> polytopeBound(clusterSize);
        while(true) {
            double unboundCounts{0.0};
            double boundCounts{0.0};
            for (size_t i = 0; i < clusterSize; ++i) {
                Transcript* transcript = transcripts[i];
                if (transcript->projectedCounts > transcript->totalCounts) {
                    transcript->projectedCounts = transcript->totalCounts;
                    polytopeBound[i] = true;
                } else if (transcript->projectedCounts < transcript->uniqueCounts) {
                    transcript->projectedCounts = transcript->uniqueCounts;
                    polytopeBound[i] = true;
                }

                if (polytopeBound[i]) {
                    boundCounts += transcript->projectedCounts;
                } else {
                    unboundCounts += transcript->projectedCounts;
                }
            }


            if (approxEqual(unboundCounts + boundCounts, clusterCounts)) {
                return;
            }

            if (unboundCounts == 0) {
                polytopeBound = boost::dynamic_bitset<>(clusterSize);
                unboundCounts = boundCounts;
                boundCounts = 0;
            }

            double normalizer = (clusterCounts - boundCounts) / unboundCounts;
            for (size_t i = 0; i < clusterSize; ++i) {
                if (!polytopeBound[i]) {
                    Transcript* transcript = transcripts[i];
                    transcript->projectedCounts *= normalizer;
                }
            }
          if (iter > maxIter) { return; }
          ++iter;
        } // end while

    }

private:
    std::list<size_t> members_;
    std::atomic<size_t> count_;
    double logMass_;
    bool active_;
};

/** A forest of transcript clusters */
class ClusterForest {
public:
    ClusterForest(size_t numTranscripts, std::vector<Transcript>& refs) :
        rank_(std::vector<size_t>(numTranscripts, 0)),
        parent_(std::vector<size_t>(numTranscripts, 0)),
        disjointSets_(&rank_[0], &parent_[0]),
        clusters_(std::vector<TranscriptCluster>(numTranscripts))
    {
        // Initially make a unique set for each transcript
        for(size_t tnum = 0; tnum < numTranscripts; ++tnum) {
            disjointSets_.make_set(tnum);
            clusters_[tnum].members_.push_front(tnum);
            clusters_[tnum].addMass(refs[tnum].mass());
        }
    }

    template <typename FragT>
    void mergeClusters(typename AlignmentBatch<FragT>::iterator start,
                       typename AlignmentBatch<FragT>::iterator finish) {
        // Use a lock_guard to ensure this is a locked (and exception-safe) operation
        std::lock_guard<std::mutex> lock(clusterMutex_);
        size_t firstCluster, otherCluster;
        auto firstTranscriptID = start->transcriptID();
        ++start;

        for (auto it = start; it != finish; ++it) {
            firstCluster = disjointSets_.find_set(firstTranscriptID);
            otherCluster = disjointSets_.find_set(it->transcriptID());
            if (otherCluster != firstCluster)  {
                disjointSets_.link(firstCluster, otherCluster);
                auto parentClust = disjointSets_.find_set(it->transcriptID());
                auto childClust = (parentClust == firstCluster)  ? otherCluster : firstCluster;
                if (parentClust == firstCluster or parentClust == otherCluster) {
                    clusters_[parentClust].merge(clusters_[childClust]);
                    clusters_[childClust].deactivate();
                } else { std::cerr << "DANGER\n"; }
            }
        }
    }


    /*
    void mergeClusters(AlignmentBatch<ReadPair>::iterator start, AlignmentBatch<ReadPair>::iterator finish) {
        // Use a lock_guard to ensure this is a locked (and exception-safe) operation
        std::lock_guard<std::mutex> lock(clusterMutex_);
        size_t firstCluster, otherCluster;
        auto firstTranscriptID = start->read1->core.tid;
        ++start;

        for (auto it = start; it != finish; ++it) {
            firstCluster = disjointSets_.find_set(firstTranscriptID);
            otherCluster = disjointSets_.find_set(it->read1->core.tid);
            if (otherCluster != firstCluster)  {
                disjointSets_.link(firstCluster, otherCluster);
                auto parentClust = disjointSets_.find_set(it->read1->core.tid);
                auto childClust = (parentClust == firstCluster)  ? otherCluster : firstCluster;
                if (parentClust == firstCluster or parentClust == otherCluster) {
                    clusters_[parentClust].merge(clusters_[childClust]);
                    clusters_[childClust].deactivate();
                } else { std::cerr << "DANGER\n"; }
            }
        }
    }
*/
    void updateCluster(size_t memberTranscript, size_t newCount, double logNewMass) {
        // Use a lock_guard to ensure this is a locked (and exception-safe) operation
        std::lock_guard<std::mutex> lock(clusterMutex_);
        auto clusterID = disjointSets_.find_set(memberTranscript);
        auto& cluster = clusters_[clusterID];
        cluster.incrementCount(newCount);
        cluster.addMass(logNewMass);
    }

    std::vector<TranscriptCluster*> getClusters() {
        std::vector<TranscriptCluster*> clusters;
        std::unordered_set<size_t> observedReps;
        for (size_t i = 0; i < clusters_.size(); ++i) {
            auto rep = disjointSets_.find_set(i);
            if (observedReps.find(rep) == observedReps.end()) {
                if (!clusters_[rep].isActive()) {
                    std::cerr << "returning a non-active cluster!\n";
                    std::exit(1);
                }
                clusters.push_back(&clusters_[rep]);
                observedReps.insert(rep);
            }
        }
        return clusters;
    }
private:
    std::vector<size_t> rank_;
    std::vector<size_t> parent_;
    boost::disjoint_sets<size_t*, size_t*> disjointSets_;
    std::vector<TranscriptCluster> clusters_;
    std::mutex clusterMutex_;
};

template <typename FragT>
void processMiniBatch(MiniBatchQueue<AlignmentGroup<FragT>>& workQueue,
                      std::condition_variable& workAvailable,
                      std::mutex& cvmutex,
                     std::vector<Transcript>& refs,
                      ClusterForest& clusterForest, std::atomic<bool>& doneParsing, std::atomic<size_t>& activeBatches,
                      tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue,
                      tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*>& alignmentGroupQueue,
                      ErrorModel& errMod,
                      FragmentLengthDistribution& fragLengthDist,
                      bool& burnedIn,
                      std::atomic<size_t>& processedReads) {

    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

    using sailfish::math::LOG_0;
    using sailfish::math::logAdd;
    using sailfish::math::logSub;

    std::chrono::microseconds sleepTime(1);
    MiniBatchInfo<AlignmentGroup<FragT>>* miniBatch;

    auto calcQuality = [&refs](bam1_t* b) -> double {
        /*
        Transcript& ref = refs[b->core.tid];
        uint32_t offset = b->core.pos;

        uint8_t* readSeq = bam1_seq(b);
        uint8_t* refSeq = (b->core.strand) ? ref.revSeq : ref.fwdSeq;

        double log10e = 1.0 / 0.43429448190325176;
        uint32_t* cigarStr = bam1_cigar(b);
        auto qualStr = bam1_qual(b);
        double qual = LOG_1;
        for(size_t i = 0; i < b->core.n_cigar; ++i) {
            switch (cigarStr[i] & 0xF) {
            case 0:
                break;
            case 8:
                qual += (static_cast<double>(-qualStr[i]) / 10.0) * log10e;
                std::cerr << "used quality for case 8\n";
                break;
            case 1:
            case 2:
                qual += std::log(cigarStr[i] >> 4) + ((static_cast<double>(-qualStr[i]) / 10.0) * log10e);
                std::cerr << "used quality for case 2\n";
                break;
            default:
                std::cerr << cigarStr[i] << "\n";
                std::exit(1);
                break;
            }
        }
        char tag[] = "NM";
        uint8_t* editDist = bam_aux_get(b, tag);
        int numMismatch{0};
        if (editDist) { numMismatch += bam_aux2i(editDist); }
        return std::exp(-numMismatch / 8.0);
        */
        return sailfish::math::LOG_1;
    };

    size_t numTranscripts = refs.size();

    while (!doneParsing) {
        miniBatch = nullptr;
        {
            std::unique_lock<std::mutex> l(cvmutex);
            workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() { return workQueue.try_pop(miniBatch) or doneParsing; });
        }
        if (miniBatch != nullptr) {
        //while(workQueue.try_pop(miniBatch)) {
            ++activeBatches;
            size_t batchReads{0};
            double logForgettingMass = miniBatch->logForgettingMass;
            std::vector<AlignmentGroup<FragT>*>& alignmentGroups = *(miniBatch->alignments);

            using TranscriptID = size_t;
            using HitIDVector = std::vector<size_t>;
            using HitProbVector = std::vector<double>;

            std::unordered_map<TranscriptID, std::vector<FragT*>> hitList;
            for (auto alnGroup : alignmentGroups) {
                for (auto& a : alnGroup->alignments()) {
                    auto transcriptID = a.transcriptID();
                    if (transcriptID < 0 or transcriptID >= refs.size()) {
                        std::cerr << "Invalid Transcript ID: " << transcriptID << "\n";
                    }
                    hitList[transcriptID].emplace_back(&a);
                }
            }

            {
                // E-step
                //boost::timer::auto_cpu_timer t(3, "E-step part 2 took %w sec.\n");

                // Iterate over each group of alignments (a group consists of all alignments reported
                // for a single read).  Distribute the read's mass proportionally dependent on the
                // current
                for (auto alnGroup : alignmentGroups) {
                    double sumOfAlignProbs{LOG_0};
                    // update the cluster-level properties
                    bool transcriptUnique{true};
                    auto firstTranscriptID = alnGroup->alignments().front().transcriptID();
                    std::unordered_set<size_t> observedTranscripts;
                    for (auto& aln : alnGroup->alignments()) {
                        auto transcriptID = aln.transcriptID();
                        auto& transcript = refs[transcriptID];
                        transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                        double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;

                        double logFragProb = sailfish::math::LOG_1;
                        //double fragLength = aln.fragLen();

                        switch (aln.fragType()) {
                            case ReadType::SINGLE_END:
                                if (aln.isLeft() and transcript.RefLength - aln.left() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(transcript.RefLength - aln.left());
                                } else if (aln.isRight() and aln.right() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(aln.right());
                                }
                                break;
                            case ReadType::PAIRED_END:
                                logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln.fragLen()));
                        }

                        // @TODO: handle this case better
                        //double fragProb = cdf(fragLengthDist, fragLength + 0.5) - cdf(fragLengthDist, fragLength - 0.5);
                        //fragProb = std::max(fragProb, 1e-3);
                        //fragProb /= cdf(fragLengthDist, refLength);

                        // Some aligners (e.g. Bowtie) don't output alignment qualities, so prepare for this.
                        /*
                        auto q1 = aln.read1->core.qual;
                        auto q2 = aln.read2->core.qual;
                        double logP1 = (q1 == 255) ? calcQuality(aln.read1) : std::log(std::pow(10.0, -q1 * 0.1));
                        double logP2 = (q2 == 255) ? calcQuality(aln.read2) : std::log(std::pow(10.0, -q2 * 0.1));
                        */

                        // The alignment probability is the product of a transcript-level term (based on abundance and) an alignment-level
                        // term below which is P(Q_1) * P(Q_2) * P(F | T)
                        double logRefLength = std::log(refLength);
                        double adjustedTranscriptLength = std::max(refLength - aln.fragLen() + 1, 1.0);
                        double logStartPosProb = std::log(1.0 / adjustedTranscriptLength);
                        // P(Fn | Tn) = Probability of selecting a fragment of this length, given the transcript is t
                        // d(Fn) / sum_x = 1^{lt} d(x)
                        //double logConditionalFragLengthProb = logFragProb  - fragLengthDist.cmf(refLength);
                        //double logProbStartPos = logStartPosProb + logConditionalFragLengthProb;
                        //double qualProb = logProbStartPos + aln.logQualProb();
                        double qualProb = -logRefLength + logFragProb + aln.logQualProb();
                        double transcriptLogCount = transcript.mass();

                        if ( transcriptLogCount != LOG_0 ) {

                            double errLike = sailfish::math::LOG_1;
                            if (burnedIn) {
                                //errLike = errMod.logLikelihood(aln, transcript);
                            }

                            //aln.logProb = (transcriptLogCount - logRefLength) + qualProb + errLike;
                            aln.logProb = transcriptLogCount + qualProb;

                            sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);
                            if (observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                                refs[transcriptID].addTotalCount(1);
                                observedTranscripts.insert(transcriptID);
                            }
                        } else {
                            aln.logProb = LOG_0;
                        }
                    }
                    // normalize the hits
                    if (sumOfAlignProbs == LOG_0) { std::cerr << "0 probability fragment; skipping\n"; continue; }
                    for (auto& aln : alnGroup->alignments()) {
                        aln.logProb -= sumOfAlignProbs;
                        auto transcriptID = aln.transcriptID();//aln.read1->core.tid;
                        auto& transcript = refs[transcriptID];
                        double r = uni(eng);
                        if (!burnedIn and r < std::exp(aln.logProb)) {
                            //errMod.update(aln, transcript, aln.logProb, logForgettingMass);
                            if (aln.fragType() == ReadType::PAIRED_END) {
                                double fragLength = aln.fragLen();//std::abs(aln.read1->core.pos - aln.read2->core.pos) + aln.read2->core.l_qseq;
                                fragLengthDist.addVal(fragLength, logForgettingMass);
                            }
                        }
                    }
                    // update the single target transcript
                    if (transcriptUnique) {
                        refs[firstTranscriptID].addUniqueCount(1);
                        clusterForest.updateCluster(firstTranscriptID, 1, logForgettingMass);
                    } else { // or the appropriate clusters
                        clusterForest.mergeClusters<FragT>(alnGroup->alignments().begin(), alnGroup->alignments().end());
                        clusterForest.updateCluster(alnGroup->alignments().front().transcriptID(), 1, logForgettingMass);
                    }

                    ++batchReads;
                } // end read group
            }// end timer

            double individualTotal = LOG_0;
            {
                // M-step
                for (auto kv = hitList.begin(); kv != hitList.end(); ++kv) {
                    auto transcriptID = kv->first;
                    // The target must be a valid transcript
                    if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr << "index " << transcriptID << " out of bounds\n"; }

                    auto& transcript = refs[transcriptID];

                    // The prior probability
                    double hitMass{LOG_0};

                    // The set of alignments that match transcriptID
                    auto& hits = kv->second;
                    std::for_each(hits.begin(), hits.end(), [&](FragT* aln) -> void {
                            if (!std::isfinite(aln->logProb)) { std::cerr << "hitMass = " << aln->logProb << "\n"; }
                            hitMass = logAdd(hitMass, aln->logProb);
                    });

                    double updateMass = logForgettingMass + hitMass;
                    individualTotal = logAdd(individualTotal, updateMass);

                    // Lock the target
                    transcript.addMass(updateMass);

                   // unlock the target
                } // end for
            } // end timer
            /*
            if (std::abs(individualTotal - clustTotal) > 1e-5) {
                std::cerr << "Cluster total = " << clustTotal << ", individual total = " << individualTotal << "\n";
            }
            */
            miniBatch->release(alignmentStructureQueue, alignmentGroupQueue);
            delete miniBatch;
            --activeBatches;
            processedReads += batchReads;
            if (processedReads >= 5000000 and !burnedIn) { burnedIn = true; }
        }
        //std::this_thread::yield();
    } // nothing left to process
}

/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatString(std::string& fmt) {
    using std::vector;
    using std::string;
    using std::map;
    using std::stringstream;

    // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
    // first convert the string to upper-case
    for (auto& c : fmt) { c = std::toupper(c); }
    // split on the delimiter ':', and put the key, value (k=v) pairs into a map
    stringstream ss(fmt);
    string item;
    map<string, string> kvmap;
    while (std::getline(ss, item, ':')) {
        auto splitPos = item.find('=', 0);
        string key{item.substr(0, splitPos)};
        string value{item.substr(splitPos+1)};
        kvmap[key] = value;
    }

    map<string, ReadType> readType = {{"SE", ReadType::SINGLE_END}, {"PE", ReadType::PAIRED_END}};
    map<string, ReadOrientation> orientationType = {{">>", ReadOrientation::SAME},
                                           {"<>", ReadOrientation::AWAY},
                                           {"><", ReadOrientation::TOWARD},
                                           {"*", ReadOrientation::NONE}};
    map<string, ReadStrandedness> strandType = {{"SA", ReadStrandedness::SA},
                                    {"AS", ReadStrandedness::AS},
                                    {"A", ReadStrandedness::A},
                                    {"S", ReadStrandedness::S},
                                    {"U", ReadStrandedness::U}};
    auto it = kvmap.find("T");
    string typeStr = "";
    if (it != kvmap.end()) {
        typeStr = it->second;
    } else {
        it = kvmap.find("TYPE");
        if (it != kvmap.end()) {
            typeStr = it->second;
        }
    }

    if (typeStr != "SE" and typeStr != "PE") {
        string e = typeStr + " is not a valid read type; must be one of {SE, PE}";
        throw std::invalid_argument(e);
    }

    ReadType type = (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
    ReadOrientation orientation = (type == ReadType::SINGLE_END) ? ReadOrientation::NONE : ReadOrientation::TOWARD;
    ReadStrandedness strandedness{ReadStrandedness::U};
    // Construct the LibraryFormat class from the key, value map
    for (auto& kv : kvmap) {
        auto& k = kv.first; auto& v = kv.second;
        if (k == "O" or k == "ORIENTATION") {
            auto it = orientationType.find(v);
            if (it != orientationType.end()) { orientation = orientationType[it->first]; } else {
                string e = v + " is not a valid orientation type; must be one of {>>, <>, ><}";
                throw std::invalid_argument(e);
            }

        }
        if (k == "S" or k == "STRAND") {
            auto it = strandType.find(v);
            if (it != strandType.end()) { strandedness = strandType[it->first]; } else {
                string e = v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
                throw std::invalid_argument(e);
            }
        }

    }
    LibraryFormat lf(type, orientation, strandedness);
    return lf;
}

/**
  *  Quantify the targets given in the file `transcriptFile` using the
  *  alignments given in the file `alignmentFile`, and write the results
  *  to the file `outputFile`.  The reads are assumed to be in the format
  *  specified by `libFmt`.
  *
  */
template <typename FragT>
void quantifyLibrary(LibraryFormat libFmt,
            bfs::path& alignmentFile, bfs::path& transcriptFile,
            bfs::path& outputFile) {

    std::string fname = alignmentFile.string();
    BAMQueue<FragT> bq(fname, libFmt);
    bq.start();

    std::atomic<size_t> i{0};
    std::vector<Transcript> refs;

    double alpha = 0.005;
    double sumAlphaInit{LOG_0};
    bam_header_t* header = bq.header();
    for (size_t i = 0; i < header->n_targets; ++i) {
        refs.emplace_back(i, header->target_name[i], header->target_len[i], alpha);
        sumAlphaInit = logAdd(sumAlphaInit, alpha);//alphas[i]);
    }

    if (!bfs::exists(transcriptFile)) {
        std::stringstream ss;
        ss << "The provided transcript file: " << transcriptFile <<
            " does not exist!\n";
        throw std::invalid_argument(ss.str());
    }

    FASTAParser fp(transcriptFile.string());
    std::cerr << "populating targets ... ";
    fp.populateTargets(refs);
    std::cerr << "done\n";

    bool burnedIn{false};

    ErrorModel errMod(1.00);

    size_t numTranscripts = refs.size();

    ClusterForest clusterForest(numTranscripts, refs);

    std::atomic<bool> doneParsing{false};
    MiniBatchQueue<AlignmentGroup<FragT>> workQueue;
    size_t maxFragLen = 800;
    size_t meanFragLen = 200;
    size_t fragLenStd = 80;
    size_t fragLenKernelN = 4;
    double fragLenKernelP = 0.5;
    FragmentLengthDistribution fragLengthDist(1.0, maxFragLen, meanFragLen, fragLenStd, fragLenKernelN, fragLenKernelP, 1);

    std::condition_variable workAvailable;
    std::mutex cvmutex;
    size_t numWorkers{6};
    std::vector<std::thread> workers;
    std::atomic<size_t> activeBatches{0};
    std::atomic<size_t> processedReads{0};
    for (size_t i = 0; i < numWorkers; ++i) {
        workers.emplace_back(processMiniBatch<FragT>, std::ref(workQueue),
                std::ref(workAvailable), std::ref(cvmutex),
                std::ref(refs), std::ref(clusterForest),
                std::ref(doneParsing), std::ref(activeBatches),
                std::ref(bq.getAlignmentStructureQueue()),
                std::ref(bq.getAlignmentGroupQueue()),
                std::ref(errMod),
                std::ref(fragLengthDist),
                std::ref(burnedIn),
                std::ref(processedReads));
    }

    size_t miniBatchSize{1000};
    size_t batchNum{0};
    size_t numProc{0};

    double logForgettingMass{std::log(1.0)};
    double forgettingFactor{0.60};

    std::vector<AlignmentGroup<FragT>*>* alignments = new std::vector<AlignmentGroup<FragT>*>;
    alignments->reserve(miniBatchSize);
    AlignmentGroup<FragT>* ag;
    while (bq.getAlignmentGroup(ag)) {
        alignments->push_back(ag);
        if (alignments->size() >= miniBatchSize) {
            ++batchNum;
            if (batchNum > 1) {
                logForgettingMass += forgettingFactor * std::log(static_cast<double>(batchNum-1)) -
                    std::log(std::pow(static_cast<double>(batchNum), forgettingFactor) - 1);
            }

            MiniBatchInfo<AlignmentGroup<FragT>>* mbi =
                new MiniBatchInfo<AlignmentGroup<FragT>>(batchNum, alignments, logForgettingMass);
            workQueue.push(mbi);
            {
                std::unique_lock<std::mutex> l(cvmutex);
                workAvailable.notify_one();
            }
            alignments = new std::vector<AlignmentGroup<FragT>*>;
            alignments->reserve(miniBatchSize);
        }
        if (numProc % 1000000 == 0) {
            const char RESET_COLOR[] = "\x1b[0m";
            char green[] = "\x1b[30m";
            green[3] = '0' + static_cast<char>(fmt::GREEN);
            char red[] = "\x1b[30m";
            red[3] = '0' + static_cast<char>(fmt::RED);
            fmt::print(stderr, "\r\r{}processed{} {} {}reads{}", green, red, numProc, green, RESET_COLOR);
        }
        ++numProc;
    }
    std::cerr << "\n";

    doneParsing = true;
    size_t tnum{0};
    for (auto& t : workers) {
        std::cerr << "killing thread" << tnum++;
        {
            std::unique_lock<std::mutex> l(cvmutex);
            workAvailable.notify_all();
        }
        t.join(); std::cerr << " done\n";
    }

    std::cerr << "writing output \n";

    std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(outputFile.c_str(), "w"), std::fclose);

    fmt::print(output.get(), "# Salmon-Read v 0.01\n");
    fmt::print(output.get(), "# ClusterID\tName\tLength\tFPKM\tNumReads\n");

    const double logBillion = std::log(1000000000.0);
    const double logNumFragments = std::log(static_cast<double>(numProc));
    auto clusters = clusterForest.getClusters();
    size_t clusterID = 0;
    for(auto cptr : clusters) {
        /*
           std::vector<double> readCounts;
           std::vector<double> posterior;
           */

        double logClusterMass = cptr->logMass();
        double logClusterCount = std::log(static_cast<double>(cptr->numHits()));

        if (logClusterMass == LOG_0) {
            std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
        }

        bool requiresProjection{false};

        auto& members = cptr->members();
        size_t clusterSize{0};
        for (auto transcriptID : members) {
            Transcript& t = refs[transcriptID];
            t.uniqueCounts = t.uniqueCount();
            t.totalCounts = t.totalCount();
            //clusterCount += t.totalCounts;
        }

        for (auto transcriptID : members) {
            Transcript& t = refs[transcriptID];
            double logTranscriptMass = t.mass();
            double logClusterFraction = logTranscriptMass - logClusterMass;
            /*
               posterior.push_back(clusterFraction);
               readCounts.push_back(clusterCount * clusterFraction);
               */
            t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
            requiresProjection |= t.projectedCounts > static_cast<double>(t.totalCounts) or
                t.projectedCounts < static_cast<double>(t.uniqueCounts);
            ++clusterSize;
        }

        if (clusterSize > 1 and requiresProjection) {
            cptr->projectToPolytope(refs);
        }

        // Now posterior has the transcript fraction
        size_t idx = 0;
        for (auto transcriptID : members) {
            auto& transcript = refs[transcriptID];
            double logLength = std::log(transcript.RefLength);
            double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
            double count = refs[transcriptID].projectedCounts;
            double countTotal = refs[transcriptID].totalCounts;
            double countUnique = refs[transcriptID].uniqueCounts;
            double fpkm = count > 0 ? fpkmFactor * count : 0.0;
            fmt::print(output.get(),
                       "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                       clusterID, transcript.RefName, transcript.RefLength,
                       fpkm, countTotal, countUnique, count, transcript.mass());

            ++idx;
        }

        ++clusterID;
    }

    fmt::print(stdout, "{}\n", fragLengthDist.toString());
}

int main(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    po::options_description generic("Sailfish read-quant options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("libtype,l", po::value<std::string>(), "Format string describing the library type")
    ("alignments,a", po::value<std::string>(), "input alignment (BAM) file")
    ("targets,t", po::value<std::string>(), "FASTA format file containing target transcripts")
    ("output,o", po::value<std::string>(), "Output quantification file");

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(generic).run();

        po::store(orderedOptions, vm);

        if (vm.count("help")) {
            std::cout << "Sailfish read-quant\n";
            std::cout << generic << std::endl;
            std::exit(0);
        }
        po::notify(vm);

        for (auto& opt : orderedOptions.options) {
            std::cerr << "[ " << opt.string_key << "] => {";
            for (auto& val : opt.value) {
                std::cerr << " " << val;
            }
            std::cerr << " }\n";
        }

        bfs::path alignmentFile(vm["alignments"].as<std::string>());
        std::cerr << vm["alignments"].as<std::string>() << "!!!!\n";
        if (!bfs::exists(alignmentFile)) {
            std::stringstream ss;
            ss << "The provided alignment file: " << alignmentFile <<
                " does not exist!\n";
            throw std::invalid_argument(ss.str());
        }

        std::string libFmtStr = vm["libtype"].as<std::string>();
        LibraryFormat libFmt = parseLibraryFormatString(libFmtStr);
        if (libFmt.check()) {
            std::cerr << libFmt << "\n";
        } else {
            std::stringstream ss;
            ss << libFmt << " is invalid!";
            throw std::invalid_argument(ss.str());
        }

        std::string fname(alignmentFile.string());
        bfs::path outputFile(vm["output"].as<std::string>());
        bfs::path transcriptFile(vm["targets"].as<std::string>());

        switch (libFmt.type) {
            case ReadType::SINGLE_END:
                quantifyLibrary<UnpairedRead>(libFmt,
                        alignmentFile, transcriptFile, outputFile);
                break;
            case ReadType::PAIRED_END:
                quantifyLibrary<ReadPair>(libFmt,
                        alignmentFile, transcriptFile, outputFile);
                break;
            default:
                std::cerr << "Cannot quantify library of unknown " <<
                        "format " << libFmt << "\n";
                std::exit(1);
        }

    } catch (po::error& e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "============\n";
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << "============\n";
        std::cerr << argv[0] << " read-quant was invoked improperly.\n";
        std::cerr << "For usage information, " <<
            "try " << argv[0] << " read-quant --help\nExiting.\n";
    }
    return 0;
}

/*
  AlignmentGroup* nag;
  bq.getAlignmentGroup(nag);

  constexpr char samToChar[] = {'*', 'A', 'C', '*', 'G', '*', '*', '*', 'T',
  '*', '*', '*', '*', '*', '*', '*'};

  for(auto& aln : nag->alignments()) {
  auto transcriptID = aln.read1->core.tid;
  auto& ref = refs[transcriptID];
  std::cerr << "ref = " << ref.RefName << "\n";

  std::cerr << "read1 = " << bam1_qname(aln.read1) << ", p = " << aln.read1->core.pos << "\n";
  size_t linit = aln.read1->core.pos;
  for (size_t offset = 0; offset < aln.read1->core.l_qseq; ++offset) {
  //std::cerr << linit + offset << " ";
  std::cerr << samToChar[ref.baseAt(linit + offset)];
  }
  std::cerr << "\n";

  uint8_t* seq = bam1_seq(aln.read1);
  for (size_t offset = 0; offset < aln.read1->core.l_qseq; ++offset) {
  //std::cerr << linit + offset << " ";
  std::cerr << samToChar[bam1_seqi(seq, offset)];
  }
  std::cerr << "\n";

  std::cerr << "read2 = " << bam1_qname(aln.read2) << ", p = " << aln.read2->core.pos << "\n";
  size_t rinit = aln.read2->core.pos;
  for (size_t offset = 0; offset < aln.read2->core.l_qseq; ++offset) {
  //std::cerr << linit + offset << " ";
  std::cerr << samToChar[ref.baseAt(rinit + offset)];
  }
  std::cerr << "\n";

  seq = bam1_seq(aln.read2);
  for (size_t offset = 0; offset < aln.read2->core.l_qseq; ++offset) {
  //std::cerr << linit + offset << " ";
  std::cerr << samToChar[bam1_seqi(seq, offset)];
  }
  std::cerr << "\n";


  }
  std::exit(1);
*/

