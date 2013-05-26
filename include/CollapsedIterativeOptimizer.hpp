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


#ifndef COLLAPSED_ITERATIVE_OPTIMIZER_HPP
#define COLLAPSED_ITERATIVE_OPTIMIZER_HPP

#include <algorithm>
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

/** Boost Includes */
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/framework/accumulator_set.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>

//#include <Eigen/Core>
#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"


#include "BiasIndex.hpp"
#include "ezETAProgressBar.hpp"
#include "LookUpTableUtils.hpp"
//#include "poisson_solver.hpp"
//#include "matrix_tools.hpp"
//#include "CompatibilityGraph.hpp"
//#include "nnls.h"

template <typename ReadHash>
class CollapsedIterativeOptimizer {

private:
    /**
    * Typedefs
    */
    typedef uint32_t TranscriptID;
    typedef uint64_t KmerID;
    typedef double KmerQuantity;
    typedef double Promiscutity;
    typedef tbb::concurrent_unordered_map< uint64_t, tbb::concurrent_vector<uint32_t> > KmerMap;

    struct TranscriptGeneVectors;
    typedef std::vector<TranscriptID> TranscriptIDVector;
    typedef std::vector<TranscriptIDVector> KmerIDMap;


    typedef std::tuple<TranscriptID, std::vector<KmerID>> TranscriptKmerSet;
    typedef std::string *StringPtr;
    typedef uint64_t TranscriptScore;
    typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> HashArray;
    typedef size_t ReadLength;

    // Necessary forward declaration
    struct TranscriptData;
    typedef std::tuple<TranscriptScore, TranscriptID> HeapPair;
    typedef typename boost::heap::fibonacci_heap<HeapPair>::handle_type Handle;

    struct TranscriptGeneVectors {
        tbb::concurrent_vector<uint32_t> transcripts;
        tbb::concurrent_vector<uint32_t> genes;
    };

    struct TranscriptData {
        TranscriptID id;
        StringPtr header;
        std::map<KmerID, KmerQuantity> binMers;
        KmerQuantity mean;
        size_t length;
    };

    struct TranscriptInfo {
        btree::btree_map<KmerID, KmerQuantity> binMers;
        KmerQuantity mean;
        ReadLength length;
        ReadLength effectiveLength;
    };

    // This struct represents a "job" (transcript) that needs to be processed
    struct TranscriptJob {
        StringPtr header;
        StringPtr seq;
        TranscriptID id;
    };

    struct TranscriptResult {
        TranscriptData *data;
        TranscriptKmerSet *ks;
    };

    struct BinmerUpdates {
        std::vector<KmerID> zeroedBinmers;
        std::vector<KmerID> updatedBinmers;
    };

    uint32_t numThreads_;
    size_t merLen_;
    ReadHash & readHash_;
    BiasIndex& biasIndex_;

    // The number of occurences above whcih a kmer is considered promiscuous
    size_t promiscuousKmerCutoff_ {std::numeric_limits<size_t>::max()};

    // Map each kmer to the set of transcripts it occurs in
    KmerIDMap transcriptsForKmer_;

    // The actual data for each transcript
    std::vector<TranscriptInfo> transcripts_;

    TranscriptGeneMap& transcriptGeneMap_;

    tbb::concurrent_unordered_set<uint64_t> genePromiscuousKmers_;

    std::vector<Promiscutity> kmerGroupPromiscuities_;
    std::vector<Promiscutity> kmerGroupBiases_;
    std::vector<KmerQuantity> kmerGroupCounts_;

    /**
     * Compute the "Inverse Document Frequency" (IDF) of a kmer within a set of transcripts.
     * The inverse document frequency is the log of the number of documents (i.e. transcripts)
     * divded by the number of documents containing this term (i.e. kmer).
     * @param  k [kmer id for which the IDF should be computed]
     * @return   [IDF(k)]
     */
    inline double _idf( uint64_t k ) {
        double df = transcriptsForKmer_[k].size();
        return (df > 0.0) ? std::log(transcripts_.size() / df) : 0.0;
    }

    /**
     * Returns true if this kmer should be considered in our estimation, false 
     * otherwise
     * @param  mer [kmer id to test for consideration]
     * @return     [true if we consider the kmer with this id, false otherwise]
     */
    inline bool _considered( uint64_t mer ) {
        // The kmer is only considered if it exists in the transcript set
        // (i.e. it's possible to cover) and it's less prmiscuous than the
        // cutoff.
        return true;
    }

    /**
     * The weight attributed to each appearence of the kmer with the given ID.
     * If the kmer with ID k occurs in m different transcripts, then 
     * _weight(k) = 1 / m.
     * @param  k [The ID of the kmer whose weight is to be computed]
     * @return   [The weight of each appearance of the kmer with the given ID]
     */
    KmerQuantity _weight( KmerID k ) {
        return 1.0 / (kmerGroupPromiscuities_[k] );
    }

    KmerQuantity _computeMedian( const TranscriptInfo& ti ) {

      using namespace boost::accumulators;
      typedef accumulator_set<double, stats<tag::median(with_p_square_quantile)>> Accumulator;

      Accumulator acc;
      for (auto binmer : ti.binMers) {
        acc(binmer.second);
      }
      
      return median(acc);
    }

    /**
     * Computes the sum of kmer counts within the transcript given by ti, but clamping
     * all non-zero counts to the given quantile.  For example, if quantile was 0.25, and
     * x and y represented the 1st and 3rd quantile of kmer counts, then every nonzero count c 
     * would be transformed as c = max(x, min(y,c));
     * 
     * @param  ti       [description]
     * @param  quantile [description]
     * @return          [description]
     */
    KmerQuantity _computeSumQuantile( const TranscriptInfo& ti, double quantile ) {
        using namespace boost::accumulators;
        typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;
        KmerQuantity sum = 0.0;
        
        accumulator_t accLow(quantile_probability = quantile);
        accumulator_t accHigh(quantile_probability = 1.0-quantile);
        for ( auto binmer : ti.binMers ) {
            if ( this->genePromiscuousKmers_.find(binmer.first) == this->genePromiscuousKmers_.end() ){
                accLow(binmer.second);
                accHigh(binmer.second);
            }        
        }

        auto cutLow = p_square_quantile(accLow);
        auto cutHigh = p_square_quantile(accHigh);

        for ( auto binmer : ti.binMers ) {
            if ( this->genePromiscuousKmers_.find(binmer.first) == this->genePromiscuousKmers_.end() ){
                sum += std::min( cutHigh, std::max( cutLow, binmer.second ) );
            }
        }
        return sum;
    }

    KmerQuantity _computeSum( const TranscriptInfo& ti ) {
        KmerQuantity sum = 0.0;
        for ( auto binmer : ti.binMers ) {
          sum += kmerGroupBiases_[binmer.first] * binmer.second;
        }
        return sum;
    }

    bool _discard( const TranscriptInfo& ti) {
        if ( ti.mean == 0.0 ) { 
            return false; 
        } else {
            ti.mean = 0.0;
            ti.binMers.clear();
            return true;
        }
    }

    KmerQuantity _computeMean( const TranscriptInfo& ti ) {
        return (ti.effectiveLength > 0.0) ? (_computeSum(ti) / ti.effectiveLength) : 0.0;
    }

    KmerQuantity _computeWeightedMean( const TranscriptInfo& ti ) {
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::count, tag::weighted_mean>, double> acc;

        for ( auto binmer : ti.binMers ) {
          if ( this->genePromiscuousKmers_.find(binmer.first) == this->genePromiscuousKmers_.end() ){
            acc(binmer.second, weight=kmerGroupBiases_[binmer.first] * _weight(binmer.first));
          }
        }

        auto nnz = count(acc);
        
        if ( nnz < ti.effectiveLength ) {
            acc(0.0, weight=ti.effectiveLength-nnz);
        }
       
        auto sum = sum_of_weights(acc);
        return sum > 0.0 ? weighted_mean(acc) : 0.0;
    }

    double _effectiveLength( const TranscriptInfo& ts ) {
        double length = 0.0;
        for ( auto binmer : ts.binMers ) {
            length += _weight(binmer.first);
        }
        return length;
    }

    template <typename T>
    T dotProd_(std::vector<T>& u, std::vector<T>& v) {

      auto dot = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), v.size()),
            T(0.0),  // identity element for summation
            [&]( const tbb::blocked_range<size_t>& r, T current_sum ) -> T {
             for (size_t i=r.begin(); i!=r.end(); ++i) {
               current_sum += (u[i]*v[i]);
             }
             return current_sum; // body returns updated value of the accumulator
             },
             []( double s1, double s2 ) {
                return s1+s2;       // "joins" two accumulated values
      });

      return dot;
    }

    void normalizeTranscriptMeans_(){
        //auto sumMean = 0.0;
        //for ( auto ti : transcripts_ ) { sumMean += ti.mean; }

        auto sumMean = tbb::parallel_reduce(
            tbb::blocked_range<size_t>(size_t(0), transcripts_.size()),
            double(0.0),  // identity element for summation
            [&, this]( const tbb::blocked_range<size_t>& r, double current_sum ) -> double {
                 for (size_t i=r.begin(); i!=r.end(); ++i) {
                     double x = this->transcripts_[i].mean;
                     current_sum += x;
                 }
                 return current_sum; // body returns updated value of the accumulator
             },
             []( double s1, double s2 ) {
                 return s1+s2;       // "joins" two accumulated values
             });



        // compute the new mean for each transcript
        tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
            [this, sumMean]( size_t tid ) -> void { this->transcripts_[tid].mean /= sumMean; });

    }

    template <typename T>
    T psumAccumulate(const tbb::blocked_range<T*>& r, T value) {
            return std::accumulate(r.begin(),r.end(),value);
    }

    template <typename T>
    T psum_(std::vector<T>& v) {
      auto func = std::bind( std::mem_fn(&CollapsedIterativeOptimizer<ReadHash>::psumAccumulate<T>), 
                             this, std::placeholders::_1, std::placeholders::_2 );
      auto sum = tbb::parallel_reduce(
        tbb::blocked_range<T*>(&v[0], &v[v.size()]),
          T{0},  // identity element for summation
          func,
          std::plus<T>()
        );
      return sum;
    }

    template <typename T>
    T pdiff_(std::vector<T>& v0, std::vector<T>& v1) {
        auto diff = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), v0.size()),
          double(0.0),  // identity element for difference
          [&]( const tbb::blocked_range<size_t>& r, double currentDiff ) -> double {
            for (size_t i=r.begin(); i!=r.end(); ++i) {
              currentDiff += v0[i] - v1[i];
            }
            return currentDiff; // body returns updated value of the accumulator
          },
          []( double s1, double s2 ) {
               return s1+s2;       // "joins" two accumulated values
          }
        );

        return diff;
    }

    template <typename T>
    T pabsdiff_(std::vector<T>& v0, std::vector<T>& v1) {
        auto diff = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(size_t(0), v0.size()),
          double(0.0),  // identity element for difference
          [&]( const tbb::blocked_range<size_t>& r, double currentDiff ) -> double {
            for (size_t i=r.begin(); i!=r.end(); ++i) {
              currentDiff = std::abs(v0[i] - v1[i]);
            }
            return currentDiff; // body returns updated value of the accumulator
          },
          []( double s1, double s2 ) {
               return s1+s2;       // "joins" two accumulated values
          }
        );

        return diff;
    }



    void normalize_(std::vector<double>& means) {
        auto sumMean = psum_(means);
        auto invSumMean = 1.0 / sumMean;

        // compute the new mean for each transcript
        tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
            [&means, invSumMean]( size_t tid ) -> void { means[tid] *= invSumMean; });

    }

    double averageCount( const TranscriptInfo& ts){
        if ( ts.binMers.size() == 0 ) { return 0.0; }
        double sum = 0.0;
        for ( auto binmer : ts.binMers ) {
            sum += kmerGroupBiases_[binmer.first] * binmer.second;
        }
        return sum / ts.binMers.size();

    }

    /**
     * This function should be run only after <b>after</b> an EM loop.
     * It estimates kmer specific biases based on how much a kmer's count
     * deviates from it's corresponding transcript's mean 
     */
    void computeKmerFidelities_() {

        std::vector<double> transcriptFidelities(transcripts_.size(), 0.0);
        tbb::parallel_for( size_t(0), transcripts_.size(),
            [this, &transcriptFidelities]( TranscriptID tid ) {
                double sumDiff = 0.0;
                auto ts = this->transcripts_[tid];
                for ( auto& b : ts.binMers ) {
                    auto diff = (this->kmerGroupBiases_[b.first] * b.second) - ts.mean;
                    sumDiff += diff*diff;
                }
                transcriptFidelities[tid] = std::sqrt(sumDiff / ts.binMers.size());
            });
        
        auto totalFidelity = std::accumulate(transcriptFidelities.begin(), transcriptFidelities.end(), 0.0);
        auto maxFidelity = *std::max_element(transcriptFidelities.begin(), transcriptFidelities.end());
        auto averageFidelity = totalFidelity / transcriptFidelities.size();

        std::cerr << "max fidelity = " << maxFidelity << "\n";

        tbb::parallel_for( size_t(0), transcriptsForKmer_.size(),
            [this, &transcriptFidelities, averageFidelity, maxFidelity]( KmerID kid ) -> void {
                double sumT = 0.0; double sumK = 0.0;
                double confidence = 0.0;
                for( auto tid : this->transcriptsForKmer_[kid] ) {
                    sumT += this->transcripts_[tid].mean;
                    sumK += this->transcripts_[tid].binMers[kid];   
                    confidence += transcriptFidelities[tid];// / averageFidelity;
                }

                confidence /= this->transcriptsForKmer_[kid].size();
                double alpha = 0.5; //std::min( 1.0, 10.0*averageFidelity / confidence);
                double bias = sumK > 0.0 ? (sumT / sumK) : 0.0;
                double prevBias = this->kmerGroupBiases_[kid];
                this->kmerGroupBiases_[kid] = alpha * bias + (1.0 - alpha) * prevBias;
            }
        );

    }

    double logLikelihood_(std::vector<double>& means) {

      std::vector<double> likelihoods(means.size(), 0.0);

        // Compute the log-likelihood 
        tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
          // for each transcript
          [&likelihoods, &means, this]( size_t tid ) {
            auto& ti = transcripts_[tid];

            double transcriptLikelihood = 0.0;
            double relativeAbundance = means[tid];
            // For each kmer in this transcript
            for ( auto& binmer : ti.binMers ) {
              likelihoods[tid] += relativeAbundance * std::log(binmer.second / this->kmerGroupCounts_[binmer.first]);
            }

          });

      return psum_(likelihoods);

    }

    // Since there's no built in hash code for vectors
  template <typename T>
  class my_hasher{
  public:
    size_t operator()(const T& x) const {
        if (x.size() == 0 ) { return 0; }
        size_t seed = x[0];
        for (auto i : boost::irange(size_t(1), x.size())) {
            boost::hash_combine(seed, static_cast<size_t>(x[i]));
        }
        return seed;
    }
  };

/**
 * Collapses all kmer which share the same transcript multiset.  Such kmers can be
 * treated as a "batch" with a count whose value is the sum of individual "batch"
 * members.
 * @param  isActiveKmer       [A bitvector which designates, for each kmer,
 *                             whether or not that kmer is active in the current
 *                             read set.]
 */
 void collapseKmers_( boost::dynamic_bitset<>& isActiveKmer ) {

    auto numTranscripts = transcriptGeneMap_.numTranscripts();

    /**
     * Map from a vector of transcript IDs to the list of kmers that have this
     * transcript list.  This allows us to collapse all kmers that exist in the
     * exact same set of transcripts into a single kmer group.
     */
    tbb::concurrent_unordered_map< TranscriptIDVector, 
                                   tbb::concurrent_vector<KmerID>,
                                   my_hasher<std::vector<TranscriptID>> > m;

     // Asynchronously print out the progress of our hashing procedure                              
     std::atomic<size_t> prog{0};
     std::thread t([this, &prog]() {
        ez::ezETAProgressBar pb(this->transcriptsForKmer_.size());
        pb.start();
        size_t prevProg{0};
        while ( prevProg < this->transcriptsForKmer_.size() ) {
            if (prog > prevProg) {
                auto diff = prog - prevProg;
                pb += diff;
                prevProg += diff;
            }
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        if (!pb.isDone()) { pb.done(); }
     });

     //For every kmer, compute it's kmer group.
     tbb::parallel_for( size_t(0), transcriptsForKmer_.size(),
        [&]( size_t j ) {
          if (isActiveKmer[j]) {
            m[ transcriptsForKmer_[j]  ].push_back(j);
        }
        ++prog;
     });

     // wait for the parallel hashing to finish
     t.join();

     std::cerr << "Out of " << transcriptsForKmer_.size() << " potential kmers, "
               << "there were " << m.size() << " distinct groups\n";

     size_t totalKmers = 0;
     size_t index = 0;
     std::vector<KmerQuantity> kmerGroupCounts(m.size());
     std::vector<Promiscutity> kmerGroupPromiscuities(m.size());
     std::vector<TranscriptIDVector> transcriptsForKmer(m.size());

     using namespace boost::accumulators;
     std::cerr << "building collapsed transcript map\n";
     for ( auto& kv : m ) {

        // For each transcript covered by this kmer group, add this group to the set of kmer groups contained in 
        // the transcript.  For efficiency, we also compute the kmer promiscuity values for each kmer
        // group here --- the promiscuity of a kmer group is simply the number of distinct transcripts in
        // which this group of kmers appears.
        auto prevTID = std::numeric_limits<TranscriptID>::max();
        KmerQuantity numDistinctTranscripts = 0.0;
        for ( auto& tid : kv.first ) {
          transcripts_[tid].binMers[index] += 1;
          // Since the transcript IDs are sorted we just have to check
          // if this id is different from the previous one
          if (tid != prevTID) { numDistinctTranscripts += 1.0; }
          prevTID = tid;
        }
        // Set the promiscuity and the set of transcripts for this kmer group
        kmerGroupPromiscuities[index] = numDistinctTranscripts;
        transcriptsForKmer[index] = kv.first;

        // Aggregate the counts attributable to each kmer into its repective
        // group's counts.
        for (auto kid : kv.second) {
            kmerGroupCounts[index] += readHash_.atIndex(kid);
        }

        // Update the total number of kmers we're accounting for
        // and the index of the current kmer group.
        totalKmers += kv.second.size();
        ++index;
      }

      std::cerr << "Verifying that the unique set encodes " << totalKmers << " kmers\n";
      std::cerr << "collapsedCounts.size() = " << transcriptsForKmer.size() << "\n";

      // update the relevant structures holding info for the full kmer
      // set with those holding the info for our collapsed kmer sets
      std::swap(kmerGroupPromiscuities, kmerGroupPromiscuities_);
      std::swap(kmerGroupCounts, kmerGroupCounts_);
      std::swap(transcriptsForKmer, transcriptsForKmer_);

      /*
      uint64_t groupCounts = 0;
      for(auto c : kmerGroupCounts_) { groupCounts += c; }
      auto tmp = psum_(kmerGroupCounts_);
      std::cerr << "groupCount(" << groupCounts << ") - parallelCount(" << tmp << ") = " << groupCounts - tmp << "\n";
      std::atomic<uint64_t> individualCounts{0};
      tbb::parallel_for(size_t{0}, readHash_.size(),
        [&](size_t i) { 
          if (isActiveKmer[i]) {
            individualCounts += readHash_.atIndex(i);
          } });
      auto diff = groupCounts - individualCounts;
      std::cerr << "groupTotal(" << groupCounts << ") - totalNumKmers(" << individualCounts << ") = " << diff << "\n";
      */
  }

  /**
   * This function should be called before performing any optimization procedure.
   * It builds all of the necessary data-structures which are used during the transcript
   * quantification procedure.
   * @param  klutfname [The name of the file containing the kmer lookup table.]
   * @param  tlutfname [The name of the file containing the transcript lookup table.]
   */
    void initialize_(
        const std::string& klutfname,
        const std::string& tlutfname,
        const bool discardZeroCountKmers) {

        // So we can concisely identify each transcript
        TranscriptID transcriptIndex {0};

        size_t numTranscripts = transcriptGeneMap_.numTranscripts();
        size_t numKmers = readHash_.size();
        auto merSize = readHash_.kmerLength();        

        size_t numActors = numThreads_;
        std::vector<std::thread> threads;

        transcripts_.resize(transcriptGeneMap_.numTranscripts());

        tbb::parallel_for( size_t{0}, transcriptGeneMap_.numTranscripts(),
            [this](size_t tid) -> void {
             this->transcripts_[tid] = TranscriptInfo {
                btree::btree_map<KmerID, KmerQuantity>(),
                0.0,
                0,
                0
                };
            }
        );

        // Get the kmer look-up-table from file
        LUTTools::readKmerLUT(klutfname, transcriptsForKmer_);

        boost::dynamic_bitset<> isActiveKmer(numKmers);

        for (auto kid : boost::irange(size_t{0}, numKmers)) {
          if ( !discardZeroCountKmers or readHash_.atIndex(kid) != 0) {
            isActiveKmer[kid] = 1;
          }
        }
        // DYNAMIC_BITSET is *NOT* concurrent?!?!
        // determine which kmers are active
        // tbb::parallel_for(size_t(0), numKmers, 
        //     [&](size_t kid) {  
        //       if ( !discardZeroCountKmers or readHash_.atIndex(kid) != 0) {
        //         isActiveKmer[kid] = 1;
        //       }
        // });

        // compute the equivalent kmer sets
        std::cerr << "\n";
        collapseKmers_(isActiveKmer);

        // we have no biases currently
        kmerGroupBiases_.resize(transcriptsForKmer_.size(), 1.0);

        // Get transcript lengths
        std::ifstream ifile(tlutfname, std::ios::binary);
        size_t numRecords {0};
        ifile.read(reinterpret_cast<char *>(&numRecords), sizeof(numRecords));
        std::cerr << "Transcript LUT contained " << numRecords << " records\n";
        for (auto i : boost::irange(size_t(0), numRecords)) {
            auto ti = LUTTools::readTranscriptInfo(ifile);
            // copy over the length, then we're done.
            transcripts_[ti->transcriptID].length = ti->length;
            transcripts_[ti->transcriptID].effectiveLength = ti->length - merSize + 1;
        }
        ifile.close();
        // --- done ---

       // tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
       //     [this]( size_t idx ) { 
       //         auto& transcript = this->transcripts_[idx];
       //         transcript.effectiveLength = transcript.effectiveLength - transcript.binMers.size();
       //         for (auto binmer : transcript.binMers) {
       //          transcript.effectiveLength += this->_weight(binmer.first);
       //         }
       // });



        size_t numRes = 0;
        std::cerr << "\n\nRemoving duplicates from kmer transcript lists ... ";
        tbb::parallel_for( size_t(0), size_t(transcriptsForKmer_.size()),
            [&numRes, this]( size_t idx ) { 
                auto& transcripts = this->transcriptsForKmer_[idx];
                // should already be sorted -- extra check can be removed eventually
                std::is_sorted(transcripts.begin(), transcripts.end());
                // Uniqify the transcripts
                auto it = std::unique(transcripts.begin(), transcripts.end()); 
                transcripts.resize(std::distance(transcripts.begin(), it));
                ++numRes;
         });
         
         std::cerr << "done\n";

         std::cerr << "Computing kmer group promiscuity rates\n";
         /* -- done
         kmerGroupPromiscuities_.resize(transcriptsForKmer_.size());
         tbb::parallel_for( size_t{0}, kmerGroupPromiscuities_.size(),
            [this]( KmerID kid ) -> void { this->kmerGroupPromiscuities_[kid] = this->_weight(kid); }
         );
         */
        
        tbb::parallel_for(size_t{0}, transcripts_.size(),
          [&, this](size_t tid) {
            auto& ti = this->transcripts_[tid];
            for (auto& binmer : ti.binMers) {
              if (binmer.second > promiscuousKmerCutoff_) {
                ti.effectiveLength -= 1.0;
              }
            }
        });

        /**
         * gene-promiscuous kmers can never be added to a transcript's counts, so
         * it's unfair to consider them in the transcripts effective length. 
         */
        std::for_each( genePromiscuousKmers_.begin(), genePromiscuousKmers_.end(),
            [this]( KmerID kmerId ) { 
                for ( auto tid : transcriptsForKmer_[kmerId] ) {
                    transcripts_[tid].effectiveLength -= 1.0;
                }
            });

        std::cerr << "done\n";

        //return mappedReads;
    }

    void _dumpCoverage( const std::string &cfname ) {
        typedef std::string* StringPtr;

        size_t numTrans = transcripts_.size();
        size_t numProc = 0;
        std::ofstream ofile(cfname);

        ofile << "# numtranscripts_\n";
        ofile << "# transcript_name_{1} num_kmers_{1} count_1 count_2 ... count_{num_kmers}\n";
        ofile << "# ... \n";
        ofile << "# transcript_name_{numtranscripts_} num_kmers_{numtranscripts_} count_1 count_2 ... count_{num_kmers_{numtranscripts_}}\n";

        ofile << transcripts_.size() << "\n";

        std::cerr << "Dumping coverage statistics to " << cfname << "\n";


        tbb::concurrent_queue<StringPtr> covQueue;
        // boost::lockfree::queue<StringPtr> covQueue(transcripts_.size());
                
        tbb::parallel_for( size_t{0}, transcripts_.size(),
            [this, &covQueue] (size_t index) -> void {

                const auto& td = this->transcripts_[index];
                
                std::stringstream ostream;
                ostream << this->transcriptGeneMap_.transcriptName(index) << " " << td.binMers.size();
                for ( auto bm : td.binMers ) {
                    ostream << " " << bm.second;
                }
                ostream << "\n";
                std::string* ostr = new std::string(ostream.str());
                covQueue.push(ostr);
                // for boost lockfree
                // while(!covQueue.push(ostr));
            }
        );


                        
        ez::ezETAProgressBar pb(transcripts_.size());
        pb.start();

        std::string* sptr = nullptr;
        while ( numProc < numTrans ) {
            while( covQueue.try_pop(sptr) ) {
                ofile << (*sptr);
                ++pb;
                ++numProc;
                delete sptr;
            }
        }

        ofile.close();

    }

    void EMUpdate_( const std::vector<double>& meansIn, std::vector<double>& meansOut ) {
      assert(meansIn.size() == meansOut.size());

      auto reqNumJobs = transcriptsForKmer_.size();                

      std::atomic<size_t> numJobs{0};
      std::atomic<size_t> completedJobs{0};

      // Print out our progress
      auto pbthread = std::thread( 
        [&completedJobs, reqNumJobs]() -> bool {
          auto prevNumJobs = 0;
          ez::ezETAProgressBar show_progress(reqNumJobs);
          show_progress.start();
          while ( prevNumJobs < reqNumJobs ) {
            if ( prevNumJobs < completedJobs ) {
              show_progress += completedJobs - prevNumJobs;
            }
            prevNumJobs = completedJobs.load();
            boost::this_thread::sleep_for(boost::chrono::seconds(1));
          }
          if (!show_progress.isDone()) { show_progress.done(); }
          return true;
        });

        //  E-Step : reassign the kmer group counts proportionally to each transcript
        tbb::parallel_for( size_t(0), size_t(transcriptsForKmer_.size()),
          // for each kmer group
          [&completedJobs, &meansIn, this]( size_t kid ) {
            auto kmer = kid;
            if ( this->genePromiscuousKmers_.find(kmer) == this->genePromiscuousKmers_.end() ){

            // for each transcript containing this kmer group
              auto &transcripts = this->transcriptsForKmer_[kmer];
              if ( transcripts.size() > 0 ) {

                double totalMass = 0.0;
                for ( auto tid : transcripts ) {
                  totalMass += meansIn[tid];
                }

                if ( totalMass > 0.0 ) {
                  double norm = 1.0 / totalMass;
                  for ( auto tid : transcripts ) {
                    if ( meansIn[tid] > 0.0 ) {
                      this->transcripts_[tid].binMers[kmer] =
                      meansIn[tid] * norm * kmerGroupBiases_[kmer] * this->kmerGroupCounts_[kmer];
                    }
                  }
                }

              }
            }
            ++completedJobs;
          });

          // wait for all kmer groups to be processed
          pbthread.join();

          double delta = 0.0;
          double norm = 1.0 / transcripts_.size();
          size_t discard = 0;

          // M-Step
          // compute the new mean for each transcript
          tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
            [this, &meansOut]( size_t tid ) -> void {
              auto& ts = this->transcripts_[tid];
                auto tsNorm = 1.0;//(ts.effectiveLength > 0.0) ? 1.0 / std::sqrt(ts.effectiveLength) : 1.0;
                meansOut[tid] = tsNorm * this->_computeMean( ts );
          });

          normalize_(meansOut);

    }


public:
    /**
     * Construct the solver with the read and transcript hashes
     */
    CollapsedIterativeOptimizer( ReadHash &readHash, TranscriptGeneMap& transcriptGeneMap,
                                 BiasIndex& biasIndex, uint32_t numThreads ) : 
                                 readHash_(readHash), merLen_(readHash.kmerLength()), 
                                 transcriptGeneMap_(transcriptGeneMap), biasIndex_(biasIndex),
                                 numThreads_(numThreads) {}


    KmerQuantity optimize(const std::string& klutfname,
                           const std::string& tlutfname,
                           const std::string &outputFile, 
                           size_t numIt, 
                           double minMean) {

        const bool discardZeroCountKmers = true;
        initialize_(klutfname, tlutfname, discardZeroCountKmers);

        KmerQuantity globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::vector<KmerID> kmerList( transcriptsForKmer_.size(), 0 );
        size_t idx = 0;

        tbb::task_scheduler_init tbb_init(numThreads_);

        std::cerr << "Computing initial coverage estimates ... ";


        std::vector<double> means0(transcripts_.size(), 0.0);
        std::vector<double> means1(transcripts_.size(), 0.0);
        std::vector<double> means2(transcripts_.size(), 0.0);
        std::vector<double> meansPrime(transcripts_.size(), 0.0);

        std::vector<double> r(transcripts_.size(), 0.0);
        std::vector<double> v(transcripts_.size(), 0.0);

        // Compute the initial mean for each transcript
        tbb::parallel_for( size_t(0), size_t(transcriptGeneMap_.numTranscripts()),
        [this, &means0]( size_t tid ) -> void {
            auto& transcriptData = this->transcripts_[tid];

            for ( auto & kv : transcriptData.binMers ) {
                auto kmer = kv.first;
                if ( this->genePromiscuousKmers_.find(kmer) == this->genePromiscuousKmers_.end() ){
                    // count is the number of times kmer appears in transcript (tid)
                    auto count = kv.second;
                    kv.second = count * this->kmerGroupCounts_[kmer] * this->_weight(kmer);
                }
            }
            transcriptData.mean = means0[tid] = this->_computeMean(transcriptData);
            //transcriptData.mean = means0[tid] = this->_computeWeightedMean(transcriptData);
            //transcriptData.mean = means0[tid] = 1.0 / this->transcripts_.size();//this->_computeMean(transcriptData);
            //transcriptData.mean = this->_computeWeightedMean(transcriptData);
            //transcriptData.mean = distribution(generator);
            //this->_computeWeightedMean( transcriptData );
        }
        );
        normalizeTranscriptMeans_();
        normalize_(means0);

        std::cerr << "done\n";
        size_t outerIterations = 1;
        /*
        for ( size_t iter = 0; iter < numIt; ++iter ) {
          std::cerr << "EM iteraton: " << iter << "\n";
          EMUpdate_(means0, means1);
          std::swap(means0, means1);
        }
        */
        
        /**
         * Defaults for these values taken from the R implementation of
         * [SQUAREM](http://cran.r-project.org/web/packages/SQUAREM/index.html).
         */
        double minStep0, minStep, maxStep0, maxStep, mStep, nonMonotonicity;
        minStep0 = minStep = 1.0;
        maxStep = maxStep = 1.0;
        mStep = 4.0;
        nonMonotonicity = 1.0;

        double negLogLikelihoodOld = std::numeric_limits<double>::infinity();
        double negLogLikelihoodNew = std::numeric_limits<double>::infinity();

        // Right now, the # of iterations is fixed, but termination should
        // also be based on tolerance
        for ( size_t iter = 0; iter < numIt; ++iter ) {
          std::cerr << "SQUAREM iteraton [" << iter << "]\n";

          // Theta_1 = EMUpdate(Theta_0)
          std::cerr << "1/3\n";
          EMUpdate_(means0, means1);

          if (!std::isfinite(negLogLikelihoodOld)) {
            negLogLikelihoodOld = -logLikelihood_(means0);
          }

          // Theta_2 = EMUpdate(Theta_1)
          std::cerr << "2/3\n";
          EMUpdate_(means1, means2);

          // r = Theta_1 - Theta_0
          tbb::parallel_for(size_t(0), transcripts_.size(),
            [&means0, &means1, &r](size_t tid) -> void { r[tid] = means1[tid] - means0[tid];}
          );

          // v = (Theta_2 - Theta_1) - r
          tbb::parallel_for(size_t(0), transcripts_.size(),
            [&means2, &means1, &r, &v](size_t tid) -> void { v[tid] = (means2[tid] - means1[tid]) - r[tid];}
          );

          double rNorm = std::sqrt(dotProd_(r,r));
          double vNorm = std::sqrt(dotProd_(v,v));
          double alphaS = rNorm / vNorm;

          alphaS = std::max(minStep, std::min(maxStep, alphaS));
          
          tbb::parallel_for(size_t(0), transcripts_.size(),
            [&r, &v, alphaS, &means0, &meansPrime](size_t tid) -> void { 
              meansPrime[tid] = std::max(0.0, means0[tid] + 2*alphaS*r[tid] + (alphaS*alphaS)*v[tid]);
            }
          );
          
          // Stabilization step
          if ( std::abs(alphaS - 1.0) > 0.01) {
            std::cerr << "alpha = " << alphaS << ". ";
            std::cerr << "Performing a stabilization step.\n";
            EMUpdate_(meansPrime, meansPrime);
          }

          /** Check for an error in meansPrime **/

          /** If there is **/
          if (std::isfinite(nonMonotonicity)) {
            negLogLikelihoodNew = -logLikelihood_(meansPrime);
            std::cerr << "logLikelihood = " << -negLogLikelihoodNew << ", ";
          } else {
            negLogLikelihoodNew = negLogLikelihoodOld;
          }

          if (negLogLikelihoodNew > negLogLikelihoodOld + nonMonotonicity) {
            std::swap(meansPrime, means2);
            negLogLikelihoodNew = -logLikelihood_(meansPrime);
            if (alphaS == maxStep) { maxStep = std::max(maxStep0, maxStep/mStep); }
            alphaS = 1.0;
          }
          std::cerr << "alpha = " << alphaS << ", ";

          if (alphaS == maxStep) { maxStep = mStep * maxStep; }
          if (minStep < 0 and alphaS == minStep) { minStep = mStep * minStep; }
          std::swap(meansPrime, means0);//EMUpdate_(meansPrime, means0);

          double delta = pabsdiff_(means0, means1);
          std::cerr << "delta = " << delta << "\n";

          if (!std::isnan(negLogLikelihoodNew)) {
            negLogLikelihoodOld = negLogLikelihoodNew;
          }

        }
        

        /*
        for ( size_t oiter = 0; oiter < outerIterations; ++oiter ) {
            for ( size_t iter = 0; iter < numIt; ++iter ) {

                auto reqNumJobs = transcriptsForKmer_.size();                
                std::cerr << "iteraton: " << iter << "\n";

                globalError = 0.0;
                numJobs = 0;
                completedJobs = 0;

                // Print out our progress
                auto pbthread = std::thread( 
                    [&completedJobs, reqNumJobs]() -> bool {
                        auto prevNumJobs = 0;
                        ez::ezETAProgressBar show_progress(reqNumJobs);
                        show_progress.start();
                        while ( prevNumJobs < reqNumJobs ) {
                            if ( prevNumJobs < completedJobs ) {
                                show_progress += completedJobs - prevNumJobs;
                            }
                            prevNumJobs = completedJobs.load();
                            boost::this_thread::sleep_for(boost::chrono::seconds(1));
                        }
                        return true;
                    });

                //  E-Step : reassign the kmer group counts proportionally to each transcript
                tbb::parallel_for( size_t(0), size_t(transcriptsForKmer_.size()),
                    // for each kmer group
                    [&kmerList, &completedJobs, this]( size_t kid ) {
                        auto kmer = kid;
                        if ( this->genePromiscuousKmers_.find(kmer) == this->genePromiscuousKmers_.end() ){

                            // for each transcript containing this kmer group
                            auto &transcripts = this->transcriptsForKmer_[kmer];
                            if ( transcripts.size() > 0 ) {

                                double totalMass = 0.0;
                                for ( auto tid : transcripts ) {
                                    totalMass += this->transcripts_[tid].mean;
                                }

                                if ( totalMass > 0.0 ) {
                                    double norm = 1.0 / totalMass;
                                    for ( auto tid : transcripts ) {
                                        if ( this->transcripts_[tid].mean > 0.0 ) {
                                            this->transcripts_[tid].binMers[kmer] =
                                            this->transcripts_[tid].mean * norm * kmerGroupBiases_[kmer] * this->kmerGroupCounts_[kmer];
                                        }
                                    }
                                }

                            }
                        }
                        ++completedJobs;
                    });

                // wait for all kmer groups to be processed
                pbthread.join();

                // reset the job counter
                completedJobs = 0;

                double delta = 0.0;
                double norm = 1.0 / transcripts_.size();

                std::vector<KmerQuantity> prevMeans( transcripts_.size(), 0.0 );
                tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
                    [this, &prevMeans]( size_t tid ) -> void { prevMeans[tid] = this->transcripts_[tid].mean; });

                std::cerr << "\ncomputing new means ... ";
                size_t discard = 0;

                // M-Step
                // compute the new mean for each transcript
                tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
                    [this, iter, numIt, norm, minMean, &discard]( size_t tid ) -> void {
                        auto& ts = this->transcripts_[tid];
                        auto tsNorm = 1.0;//(ts.effectiveLength > 0.0) ? 1.0 / std::sqrt(ts.effectiveLength) : 1.0;
                        //ts.mean = tsNorm * this->_computeWeightedMean( ts );
                        //ts->mean = tsNorm * this->averageCount( ts );
                        ts.mean = tsNorm * this->_computeMean( ts );
                        //ts.mean = tsNorm * this->_computeMedian( ts );
                });

                normalizeTranscriptMeans_();
                for( auto tid : boost::irange(size_t{0}, prevMeans.size()) ){
                    delta += std::abs( transcripts_[tid].mean - prevMeans[tid] );
                }

                std::cerr << "done\n";
                std::cerr << "total variation in mean = " << delta << "\n";
                std::cerr << "discarded " << discard << " transcripts in this round whose mean was below " << minMean << "\n";
            }

            std::cerr << "end of outer iteration " << oiter << " recomputing biases\n";
            // Thresholding
            tbb::parallel_for( size_t(0), size_t(transcripts_.size()),
                [this, minMean](size_t tid) -> void {
                    auto& ts = this->transcripts_[tid];
                    if (ts.mean < minMean) { 
                        ts.mean = 0.0; 
                        for (auto& kv : ts.binMers) {
                            kv.second = 0.0;
                        }
                    }
            });
            //computeKmerFidelities_();
        }
        */
        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(transcripts_.size());
        pb.start();

        std::atomic<size_t> totalNumKmers{0};
        //  Count the total number of kmers (i.e. mapped reads)
        tbb::parallel_for( size_t{0}, readHash_.size(),
          [&]( size_t kid ) { totalNumKmers += this->readHash_.atIndex(kid);} );

        //std::vector<double> fracNuc(means0.size(), 0.0);
        std::vector<double> cha(means0.size(), 0.0);
        std::vector<double> fracTran(means0.size(), 0.0);
        std::vector<double> tranExp(means0.size(), 0.0);
        std::vector<double> transcriptReads(means0.size(), 0.0);
        auto estimatedGroupTotal = psum_(kmerGroupCounts_);

        // Compute transcript fraction (\tau_i in RSEM)
        tbb::parallel_for(size_t(0), transcripts_.size(), 
          [&](size_t i) { 
            auto& ts = transcripts_[i];
            cha[i] = means0[i] * ts.effectiveLength * 1.1231075637550727;
          });
        normalize_(cha);

        // Compute transcript fraction (\tau_i in RSEM)
        tbb::parallel_for(size_t(0), transcripts_.size(), 
          [&](size_t i) { 
            auto& ts = transcripts_[i];
            fracTran[i] = (ts.effectiveLength > 0.0) ? means0[i] / ts.effectiveLength : 0.0;
            tranExp[i] = estimatedGroupTotal * cha[i];
            //tranExp[i] = _computeSum(ts);
            transcriptReads[i] = _computeSum(ts);
        });
        normalize_(fracTran);

        auto estimatedTotal = psum_(tranExp);
        double estimatedReadLength = readHash_.averageLength();
        double kmersPerRead = ((estimatedReadLength - merLen_) + 1);

        double estimatedNumReads = estimatedTotal / kmersPerRead;

        std::cerr << "averageLength = " << readHash_.averageLength() << "\n";
        std::cerr << "kmersPerRead = " << kmersPerRead << "\n";
        std::cerr << "estimatedNumReads = " << estimatedNumReads << "\n";

        auto diff = estimatedTotal - totalNumKmers;
        auto diff2 = estimatedGroupTotal - totalNumKmers;
        std::cerr << "estimatedTotal(" << estimatedTotal << ") - totalNumKmers(" <<
                     totalNumKmers << ") = " << diff << "\n";

        std::cerr << "estimatedGroupTotal(" << estimatedGroupTotal << ") - totalNumKmers(" <<
                     totalNumKmers << ") = " << diff2 << "\n";

        std::ofstream ofile( outputFile );
        size_t index = 0;
        double million = std::pow(10.0, 6);
        double billion = std::pow(10.0, 9);
        ofile << "Transcript" << '\t' << "Length" << '\t' << "Effective Length" << '\t' << 
                 "TPM" << '\t' << "RPKM" << '\t' << "Weighted Mapped Reads\n";
        for ( auto i : boost::irange(size_t{0}, transcripts_.size()) ) {
          auto& ts = transcripts_[i]; 
          // expected # of kmers coming from transcript i
          auto ci = tranExp[i] * 2;
          // expected # of reads coming from transcript i
          auto ri = ci / kmersPerRead;
          auto effectiveLength = ts.length - std::floor(estimatedReadLength) + 1;
          auto rpkm = (effectiveLength > 0 and ri > 0.0) ? 
                      (billion * (ci / (estimatedGroupTotal * ts.effectiveLength))) : 0.0;
               rpkm = (rpkm < 0.01) ? 0.0 : rpkm;
          ofile << transcriptGeneMap_.transcriptName(index) << 
                   '\t' << ts.length << '\t' <<
                   ts.effectiveLength << '\t' <<
                   fracTran[i] * million << '\t' << 
                   rpkm << '\t' <<
                   ri << "\n";
 
          ++index;
          ++pb;
        }
        ofile.close();

        auto writeCoverageInfo = false;
        if ( writeCoverageInfo ) {
            std::string cfname("transcriptCoverage.txt");
            _dumpCoverage( cfname );
        }

    }


};

#endif // ITERATIVE_OPTIMIZER_HPP
