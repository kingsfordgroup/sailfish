#ifndef ITERATIVE_OPTIMIZER_HPP
#define ITERATIVE_OPTIMIZER_HPP

#include <algorithm>
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <thread>

#include "btree_map.h"

/** Boost Includes */
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>

//#include <Eigen/Core>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/task_scheduler_init.h"

//#include "nnls.h"

#include "ezETAProgressBar.hpp"

template <typename ReadHashT, typename TranscriptHashT>
class IterativeOptimizer {

private:
    /**
    * Typedefs
    */
    // A KmerMapT is a map from a kmer (encoded as an integer) to the set
    // of transcripts where it occurs
    typedef uint32_t TranscriptIDT;
    typedef uint64_t KmerIDT;
    typedef double KmerQuantityT;
    //typedef std::unordered_map< KmerIDT, std::vector<TranscriptIDT> > KmerMapT;
    typedef tbb::concurrent_unordered_map< uint64_t, tbb::concurrent_vector<uint32_t> > KmerMapT;

    struct TranscriptGeneVectors;
    typedef std::vector<TranscriptGeneVectors> KmerIDMap;


    typedef std::tuple<TranscriptIDT, std::vector<KmerIDT>> TranscriptKmerSet;
    typedef std::string *StringPtrT;
    typedef uint64_t TranscriptScoreT;
    typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> HashArrayT;
    typedef size_t ReadLengthT;

    // Necessary forward declaration
    struct TranscriptDataT;
    typedef std::tuple<TranscriptScoreT, TranscriptIDT> HeapPair;
    typedef typename boost::heap::fibonacci_heap<HeapPair>::handle_type HandleT;

    struct TranscriptGeneVectors {
        tbb::concurrent_vector<uint32_t> transcripts;
        tbb::concurrent_vector<uint32_t> genes;
    };

    struct TranscriptDataT {
        TranscriptIDT id;
        StringPtrT header;
        std::map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        size_t length;
    };

    struct TranscriptInfo {
        btree::btree_map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        ReadLengthT length;
        ReadLengthT effectiveLength;
    };

    // This struct represents a "job" (transcript) that needs to be processed
    struct TranscriptJob {
        StringPtrT header;
        StringPtrT seq;
        TranscriptIDT id;
    };

    struct TranscriptResult {
        TranscriptDataT *data;
        TranscriptKmerSet *ks;
    };

    struct BinmerUpdates {
        std::vector<KmerIDT> zeroedBinmers;
        std::vector<KmerIDT> updatedBinmers;
    };

    size_t _merLen;
    ReadHashT &_readHash;
    TranscriptHashT &_transcriptHash;

    // The number of occurences above whcih a kmer is considered promiscuous
    size_t _promiscuousKmerCutoff {10000};

    // Map each kmer to the set of transcripts it occurs in
    //KmerMapT _transcriptsForKmer;
    KmerIDMap _transcriptsForKmer;

    // The actual data for each transcript
    std::vector<TranscriptInfo *> _transcripts;

    TranscriptGeneMap& _transcriptGeneMap;

    tbb::concurrent_unordered_set<uint64_t> _genePromiscuousKmers;

    // Should the given kmer be considered?
    inline bool _considered( uint64_t mer ) {
        // The kmer is only considered if it exists in the transcript set
        // (i.e. it's possible to cover) and it's less prmiscuous than the
        // cutoff.
        return ( _readHash[mer] > 0 and
                 _transcriptHash[mer] > 0 and
                 _transcriptHash[mer] < _promiscuousKmerCutoff );
    }

    KmerQuantityT _weight( KmerIDT k ) {
        //return 1.0 / (_transcriptHash[k]);
        return 1.0 / (_transcriptHash.atIndex(k));
    }

    KmerQuantityT _computeMedian( TranscriptInfo* ti ) {

        KmerQuantityT median = 0.0;
        auto& binMers = ti->binMers;
        auto len = binMers.size();
        if ( len > 0 ) {
            if ( len % 2 == 0 ) {
                auto it = binMers.begin();
                for ( size_t i = 0; i < (len / 2); ++i, ++it ) {}
                median = it->second; ++it;
                median += it->second;
                median /= 2.0;
            } else {
                auto it = binMers.begin();
                for ( size_t i = 0; i < (len / 2) + 1; ++i, ++it ) {}
                median = it->second;
            }
        }
        return median;
    }

    KmerQuantityT _computeSum( TranscriptInfo* ti ) {
        KmerQuantityT sum = 0.0;
        for ( auto binmer : ti->binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
                sum += binmer.second;
            }
        }
        return sum;
    }

    KmerQuantityT _computeMean( TranscriptInfo* ti ) {
        return 1.0 / ti->length * _computeSum(ti);
    }

    KmerQuantityT _computeWeightedMean( TranscriptInfo* ti ) {
        double denom = 0.0;
        double sum = 0.0;
        for ( auto binmer : ti->binMers ) {
          if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
            auto w = _weight(binmer.first);
            sum += w * binmer.second;
            denom += w;
          }
        }

        return (denom > 0.0) ? (sum / denom) : 0.0;
    }

    double _effectiveLength( TranscriptInfo* ts ) {
        double length = 0.0;
        for ( auto binmer : ts->binMers ) {
            length += _weight(binmer.first);
        }
        return length;
    }

    size_t _initialize( const std::vector<std::string> &transcriptFiles ) { 

        char** fnames = new char*[transcriptFiles.size()];
        size_t z{0};
        size_t numFnames{0};
        for ( auto& s : transcriptFiles ){
            // Ugly, yes?  But this is not as ugly as the alternatives.
            // The char*'s contained in fname are owned by the transcriptFiles
            // vector and need not be manually freed.
            fnames[numFnames] = const_cast<char*>(s.c_str());
            ++numFnames;
        }

        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);

        // So we can concisely identify each transcript
        TranscriptIDT transcriptIndex {0};

        size_t numActors = 12;
        std::vector<std::thread> threads;

        _transcripts.resize( _transcriptGeneMap.numTranscripts(), nullptr );

        _transcriptsForKmer.resize( _transcriptHash.size() );

        bool done {false};
        std::atomic<size_t> numRes {0};

        // Start the thread that will print the progress bar
        threads.push_back( std::thread( [&numRes, this] () {
            auto numJobs = this->_transcriptGeneMap.numTranscripts();
            ez::ezETAProgressBar show_progress(numJobs);
            //boost::progress_display show_progress( numJobs );
            size_t lastCount = numRes;
            show_progress.start();
            while ( lastCount < numJobs ) {
                auto diff = numRes - lastCount;
                if ( diff > 0 ) {
                    show_progress += static_cast<unsigned int>(diff);
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(1010));
            }
        }) );

        std::cerr << "Processing transcripts\n";

        // Start the desired number of threads to parse the transcripts
        // and build our data structure.
        for (size_t i = 0; i < numActors - 1; ++i) {

            threads.push_back( std::thread( 
                [&numRes, &parser, this]() -> void {

                    // Each thread gets it's own stream
                    jellyfish::parse_read::read_t* read;
                    jellyfish::parse_read::thread stream = parser.new_thread();

                    while ( (read = stream.next_read()) ) {
                        ++numRes;
                        // Strip the newlines from the read and put it into a string (newSeq)
                        std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e) - 1 );
                        auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
                        auto readLen = std::distance( seq.begin(), newEnd );
                        std::string newSeq(seq.begin(), seq.begin() + readLen);

                        // Lookup the ID of this transcript in our transcript -> gene map
                        auto transcriptIndex = this->_transcriptGeneMap.findTranscriptID( 
                            std::string( read->header, read->hlen ) );

                        auto geneIndex = this->_transcriptGeneMap.gene(transcriptIndex);

                        // The set of kmers in this transcript
                        //auto ts = new TranscriptKmerSet {transcriptIndex, {}};

                        // We know how many kmers there will be, so allocate the space
                        size_t numKmers {readLen - this->_merLen + 1};
                        //std::get<1>(*ts).reserve(numKmers);
                        
                        auto procRead = new TranscriptInfo {
                            btree::btree_map<KmerIDT, KmerQuantityT>(),
                            0.0,
                            readLen,
                            0
                        };
                        // The binary representation of the kmers in this transcript
                        //btree::btree_map<KmerIDT, KmerQuantityT> binMers;

                        // Iterate over the kmers
                        ReadLengthT effectiveLength(0);
                        float weightedLen = 0.0;
                        size_t coverage = _merLen;
                        auto INVALID = this->_transcriptHash.INVALID;
                        for ( auto offset : boost::irange( size_t(0), numKmers ) ) { 
                            // the kmer and it's uint64_t representation
                            auto mer = newSeq.substr( offset, this->_merLen );                            
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
                            auto binMerId = this->_transcriptHash.id(binMer);
                            // Only count and track kmers which should be considered
                            //if ( this->_considered(binMer) ) {
                            if ( binMerId != INVALID and this->_readHash.atIndex(binMerId) > 0 ) {
                                weightedLen += (1.0 / this->_transcriptHash.atIndex(binMerId));
                                effectiveLength += coverage;
                                coverage = 1;
                                //binMers[binMer] += 1.0;
                                procRead->binMers[binMerId] += 1.0;
                                
                                //this->_transcriptsForKmer[binMer].push_back(transcriptIndex);
                                this->_transcriptsForKmer[binMerId].transcripts.push_back(transcriptIndex);
                                this->_transcriptsForKmer[binMerId].genes.push_back(geneIndex);
                            } else {
                                coverage += (coverage < _merLen) ? 1 : 0;
                            }
                        }
                        procRead->effectiveLength = effectiveLength;
                        /*
                        auto procRead = new TranscriptInfo {
                            binMers, // the set of kmers
                            0.0, // initial desired mean
                            readLen, // real transcript length
                            effectiveLength // effective transcript length
                        };*/
                        this->_transcripts[transcriptIndex] = procRead;
                    }
            }) );

        }

        // Wait for all of the threads to finish
        for ( auto& thread : threads ){ thread.join(); }

        numRes = 0;
        tbb::task_scheduler_init tbb_init;

        std::cerr << "Determining gene-promiscuous kmers ... ";
        
        //tbb::parallel_for_each( _transcriptsForKmer.begin(), _transcriptsForKmer.end(),
        tbb::parallel_for( size_t(0), size_t(_transcriptsForKmer.size()),
            [&numRes, this]( size_t idx ) { //KmerMapT::value_type& kt ) {

                /*
                //auto& transcripts = kt.second;
                */
                
                // Uniqify the transcripts
                auto& transcripts = this->_transcriptsForKmer[idx].transcripts;
                std::sort(transcripts.begin(), transcripts.end());
                auto it = std::unique(transcripts.begin(), transcripts.end()); 
                transcripts.resize( it - transcripts.begin() );
                
                
                auto numGenes = this->_transcriptsForKmer[idx].genes.size();
                if (numGenes > 1) {
                    auto geneId = this->_transcriptsForKmer[idx].genes[0];
                    for ( auto i : boost::irange(size_t(1), this->_transcriptsForKmer[idx].genes.size()) ) {
                      auto id = this->_transcriptsForKmer[idx].genes[i];
                      if( id != geneId) { this->_genePromiscuousKmers.insert(idx); break; }
                  }
                }
                


                // auto geneId = this->_transcriptGeneMap.gene(transcripts[0]);
                // for ( auto i : boost::irange(size_t(1), transcripts.size()) ) {
                //     auto id = this->_transcriptGeneMap.gene(transcripts[i]);
                //     if( id != geneId) { this->_genePromiscuousKmers.insert(kt.first); break; }
                // }

                /*
                std::unordered_set<uint32_t> geneSet;
                for ( auto t : transcripts ){
                    geneSet.insert(this->_transcriptGeneMap.gene(t));
                }
                if ( geneSet.size() > 1 ){ this->_genePromiscuousKmers.insert(kt.first); }
                */
                ++numRes;
            }
        );
        std::cerr << "done\n";

        std::cerr << "computing number of mapped (usable) reads\n";
        size_t mappedReads = 0;
        for ( auto kidx : boost::irange(size_t(0), _transcriptsForKmer.size()) ) {
            if ( _genePromiscuousKmers.find(kidx) == _genePromiscuousKmers.end() ) {
                mappedReads += _readHash.atIndex(kidx);
            }
        }

        std::cerr << "Computing initial coverage estimates ... ";

        tbb::parallel_for( size_t(0), size_t(_transcriptGeneMap.numTranscripts()),
        [this]( size_t tid ) {
            auto transcriptData = this->_transcripts[tid];
            for ( auto & kv : transcriptData->binMers ) {
                auto kmer = kv.first;
                if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                    auto count = kv.second;
                    //kv.second = count * this->_readHash[kmer] * this->_weight(kmer);
                    kv.second = count * this->_readHash.atIndex(kmer) * this->_weight(kmer);
                }
            }

            transcriptData->mean = this->_computeWeightedMean( transcriptData );
        }
        );


        std::cerr << "done\n";
        return mappedReads;
    }

    void _dumpCoverage( const std::string &cfname ) {

        std::ofstream ofile(cfname);

        ofile << "# num_transcripts\n";
        ofile << "# transcript_name_{1} num_kmers_{1} count_1 count_2 ... count_{num_kmers}\n";
        ofile << "# ... \n";
        ofile << "# transcript_name_{num_transcripts} num_kmers_{num_transcripts} count_1 count_2 ... count_{num_kmers_{num_transcripts}}\n";

        ofile << _transcripts.size() << "\n";
        size_t index = 0;
        for ( auto & td : _transcripts ) {

            ofile << _transcriptGeneMap.transcriptName(index) << " " << td->binMers.size();
            for ( auto bm : td->binMers ) {
                ofile << " " << bm.second;
            }
            ofile << "\n";
            ++index;
        }

        ofile.close();

    }


    std::vector<double> _estimateIsoformAbundance( const size_t geneID, const size_t mappedReads
      //std::vector<std::vector<uint64_t>> &transcriptKmers,
      //std::vector<double>& transcriptLengths,
      //size_t mappedReads,
      //tbb::concurrent_unordered_map<uint64_t, uint32_t>& genePromiscuousKmers
      ) {

        typedef std::vector<std::vector<double>> IsoOptMatrix;
        typedef std::vector<double> IsoOptVector;

        auto transcripts = _transcriptGeneMap.transcriptsForGene( geneID );
        size_t numTranscripts = transcripts.size();

        std::unordered_map<uint64_t, size_t> kmerIDs;
        std::unordered_set<uint64_t> geneKmers;
        // reindex the kmers for this gene
        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first;
                auto kmerIt = geneKmers.find(kmer);
                if ( kmerIt == geneKmers.end() ) {
                    kmerIDs[kmer] = geneKmers.size();
                    geneKmers.insert(kmer);
                }
            }            
        }

        size_t numKmers = geneKmers.size();
        if ( numKmers == 0 ) {
            return IsoOptVector();
        }

        // Create a numTranscripts x numKmers matrix to hold the "rates"
        IsoOptVector x(numTranscripts, 0.0);
        IsoOptMatrix A(numTranscripts, IsoOptVector(numKmers, 0.0) );
        IsoOptVector counts(numKmers, 0.0);
        double numReads = 0.0;

        for ( size_t i = 0; i < transcripts.size(); ++i ) {
            for ( auto kc : _transcripts[ transcripts[i] ]->binMers ) {
                auto kmer = kc.first; auto count = kc.second;
                size_t j = kmerIDs[kmer];
                A[i][j] = 1.0;
                if ( counts[j] == 0 ) {
                    auto promIt = _genePromiscuousKmers.find(kmer);
                    double normFact = (promIt == _genePromiscuousKmers.end() ) ? 1.0 : 0.0;
                    counts[j] = _readHash.atIndex(kmer) * normFact;
                    numReads += counts[j];
                }
            }
        }

        
        std::unique_ptr<IsoOptMatrix> collapsedA;
        std::unique_ptr<IsoOptVector> collapsedCounts;
        matrix_tools::collapseIntoCategories(A, counts, collapsedA, collapsedCounts);
        
        for (size_t i = 0; i < collapsedA->size(); ++i) {
            auto tlen = _transcripts[ transcripts[i] ]->length;
            for (size_t j = 0; j < (*collapsedA)[0].size(); ++j) {
                (*collapsedA)[i][j] = (*collapsedA)[i][j] / ( tlen * mappedReads / 1000000000.0 );
            //(*collapsedA)[i][j] *= numReads;
            }
        }

        //transcriptLengths[transcripts[i]] / 1000.0 * mappedReads / 1000000.0;

        solve_likelihood( *collapsedA, *collapsedCounts, x, false );

        //solve_likelihood( A, counts, x, false );
        //for ( auto& e : x ) { e *= numReads; }
        return x;
    }



public:
    /**
     * Construct the solver with the read and transcript hashes
     */
    IterativeOptimizer( ReadHashT &readHash, TranscriptHashT  &transcriptHash, TranscriptGeneMap& transcriptGeneMap ) :
        _readHash(readHash), _transcriptHash(transcriptHash), _merLen(transcriptHash.kmerLength()), 
        _transcriptGeneMap(transcriptGeneMap) {}


    KmerQuantityT optimize( const std::vector<std::string> &transcriptFiles, const std::string &outputFile, size_t numIt ) {

        _initialize(transcriptFiles);

        // Holds results that we will use to update the index map
        boost::lockfree::fifo<KmerIDT> q;
        size_t numActors = 1;
        boost::threadpool::pool tp(numActors);

        KmerQuantityT globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::mutex nnlsmutex;
        std::vector<KmerIDT> kmerList( _transcriptsForKmer.size(), 0 );
        size_t idx = 0;

        /* TranscriptsForKmer change 
        // For each kmer
        for ( auto & kmerTranscripts :  _transcriptsForKmer ) {
            auto kmer = kmerTranscripts.first;
            kmerList[idx] = kmer;
            ++idx;
        }
        */
       
        tbb::task_scheduler_init tbb_init;

        for ( size_t iter = 0; iter < numIt; ++iter ) {

            auto reqNumJobs = _transcriptsForKmer.size();

            std::cerr << "iteraton: " << iter << "\n";

            globalError = 0.0;
            numJobs = 0;
            completedJobs = 0;

            auto pbthread = std::thread( 
                [&completedJobs, reqNumJobs]() -> bool {
                    auto prevNumJobs = 0;
                    ez::ezETAProgressBar show_progress(reqNumJobs);
                    //boost::progress_display show_progress( reqNumJobs );
                    show_progress.start();
                    while ( prevNumJobs < reqNumJobs ) {
                        if ( prevNumJobs < completedJobs ) {
                            show_progress += completedJobs - prevNumJobs;
                        }
                        prevNumJobs = completedJobs.load();

                        boost::this_thread::sleep_for(boost::chrono::milliseconds(1010));
                    }
                    show_progress.done();
                    return true;
                }
            );

            tbb::parallel_for( size_t(0), size_t(_transcriptsForKmer.size()),
                
                [&kmerList, &completedJobs, this]( size_t kid ) {

                        auto kmer = kid;

                        if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                        //auto &transcripts = this->_transcriptsForKmer[kmer];
                        auto &transcripts = this->_transcriptsForKmer[kmer].transcripts;
                        if ( transcripts.size() > 1 ) {

                            double totalMass = 0.0;
                            for ( auto tid : transcripts ) {
                                totalMass += this->_transcripts[tid]->mean;
                            }

                            double norm = 1.0 / totalMass;
                            for ( auto tid : transcripts ) {
                                //this->_transcripts[tid]->binMers[kmer] = ( totalMass > 0.0 ) ?
                                //  norm * this->_readHash[kmer] * this->_transcripts[tid]->mean :
                                //  0.0;
                                this->_transcripts[tid]->binMers[kmer] = ( totalMass > 0.0 ) ?
                                  norm * this->_readHash.atIndex(kmer) * this->_transcripts[tid]->mean :
                                  0.0;

                            }

                        }
                        }
                        ++completedJobs;
                    }
            );

            pbthread.join();
            //tp.wait();

            // reset the job counter
            completedJobs = 0;

            double delta = 0.0;
            double norm = 1.0 / _transcripts.size();

            std::cerr << "\ncomputing new means ... ";
            // compute the new mean for each transcript
            tbb::parallel_for( size_t(0), size_t(_transcripts.size()),
                [this, iter, numIt, norm, &delta]( size_t tid ) -> void {
                        // this thread claimed myJobID;
                        auto ts = this->_transcripts[tid];
                        auto prevMean = ts->mean;
                        ts->mean = this->_computeWeightedMean( ts );
                        delta += std::abs( ts->mean - prevMean ) * norm;
                }
            );
            std::cerr << "done\n";
            std::cerr << "average variation in mean = " << delta << "\n";
        }


        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(_transcripts.size());
        pb.start();

        std::ofstream ofile( outputFile );
        size_t index = 0;
        ofile << "Transcript" << '\t' << "Length" << '\t' << "Effective Length" << '\t' << "Weighted Mapped Reads" << '\n';
        for ( auto ts : _transcripts ) {
            ofile << _transcriptGeneMap.transcriptName(index) << '\t' << ts->length << '\t' <<
                    ts->effectiveLength << '\t' << _computeSum(ts) << "\n";
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

    KmerQuantityT optimizePoisson( const std::vector<std::string> &transcriptFiles, const std::string &outputFile ) {

        auto mappedReads = _initialize(transcriptFiles);

        // Holds results that we will use to update the index map
        boost::lockfree::fifo<KmerIDT> q;
        size_t numActors = 1;
        boost::threadpool::pool tp(numActors);

        KmerQuantityT globalError {0.0};
        bool done {false};
        std::atomic<size_t> numJobs {0};
        std::atomic<size_t> completedJobs {0};
        std::mutex nnlsmutex;
        std::vector<KmerIDT> kmerList( _transcriptsForKmer.size(), 0 );
        size_t idx = 0;


        std::vector<double> abundances( _transcriptGeneMap.numTranscripts(), 0.0 );

        std::atomic<size_t> numGenesProcessed{0};

        auto abundanceEstimateProgress = std::thread( [&numGenesProcessed, this] () -> void {
            auto numJobs = this->_transcriptGeneMap.numGenes();
            ez::ezETAProgressBar pbar(numJobs);
            pbar.start();
            size_t lastCount = numGenesProcessed;
            while ( lastCount < numJobs ) {
                auto diff = numGenesProcessed - lastCount;
                if ( diff > 0 ) {
                    pbar += diff;
                    lastCount += diff;
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(1010));
            }
        };
        std::cerr << "Processing abundances\n";

        // process each gene and solve the estimation problem
        tbb::parallel_for( size_t {0}, _transcriptGeneMap.numGenes(), size_t {1},
            [this]( const size_t & geneID ) {

                auto abundance = this->_estimateIsoformAbundance(geneID, mappedReads);
                auto transcripts = this->_transcriptGeneMap.transcriptsForGene(geneID);

                for ( size_t i = 0; i < abundance.size(); ++i ) {
                    abundances[ transcripts[i] ] = abundance[i];
                }
                ++numGenesProcessed;
            }
        );

        abundanceEstimateProgress.join();

        std::cerr << "Writing output\n";
        ez::ezETAProgressBar pb(_transcripts.size());
        pb.start();

        std::ofstream ofile( outputFile );
        size_t index = 0;
        ofile << "Transcript" << '\t' << "Length" << '\t' << "Effective Length" << '\t' << "Abundance" << '\n';
        for ( auto ts : _transcripts ) {
            ofile << _transcriptGeneMap.transcriptName(index) << '\t' << ts->length << '\t' <<
                    ts->effectiveLength << '\t' << abundances[index] << "\n";
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
