#ifndef ITERATIVE_OPTIMIZER_HPP
#define ITERATIVE_OPTIMIZER_HPP

#include <algorithm>
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <mutex>
#include <thread>

/** Boost Includes */
#include <boost/program_options.hpp>
#include <boost/progress.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <Eigen/Core>

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

#include "nnls.h"

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
    typedef std::tuple<TranscriptIDT, std::vector<KmerIDT>> TranscriptKmerSet;
    typedef std::string *StringPtrT;
    typedef uint64_t TranscriptScoreT;
    typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> HashArrayT;

    // Necessary forward declaration
    struct TranscriptDataT;
    typedef std::tuple<TranscriptScoreT, TranscriptIDT> HeapPair;
    typedef typename boost::heap::fibonacci_heap<HeapPair>::handle_type HandleT;

    struct TranscriptDataT {
        TranscriptIDT id;
        StringPtrT header;
        std::map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        size_t length;
    };

    struct TranscriptInfo {
        std::map<KmerIDT, KmerQuantityT> binMers;
        KmerQuantityT mean;
        size_t length;
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
    KmerMapT _transcriptsForKmer;

    // The actual data for each transcript
    std::vector<TranscriptInfo *> _transcripts;

    TranscriptGeneMap& _transcriptGeneMap;

    tbb::concurrent_unordered_set<uint64_t> _genePromiscuousKmers;

    // Should the given kmer be considered?t
    inline bool _considered( uint64_t mer ) {
        // The kmer is only considered if it exists in the transcript set
        // (i.e. it's possible to cover) and it's less prmiscuous than the
        // cutoff.
        return ( _readHash[mer] > 0 and
                 _transcriptHash[mer] > 0 and
                 _transcriptHash[mer] < _promiscuousKmerCutoff );
    }

    KmerQuantityT _weight( KmerIDT k ) {
        return 1.0 / (1.0 + _transcriptHash[k]);
    }

    KmerQuantityT _computeMean( std::map<KmerIDT, KmerQuantityT> &binMers ) {
        KmerQuantityT mean = 0.0;
        auto norm = 1.0 / binMers.size();
        for ( auto binmer : binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
            mean += norm * binmer.second;
            }
        }
        return mean;
    }

    KmerQuantityT _computeMedian( std::map<KmerIDT, KmerQuantityT> &binMers ) {
        KmerQuantityT median = 0.0;
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

    KmerQuantityT _computeSum( std::map<KmerIDT, KmerQuantityT> &binMers ) {
        KmerQuantityT mean = 0.0;
        for ( auto binmer : binMers ) {
            if ( this->_genePromiscuousKmers.find(binmer.first) == this->_genePromiscuousKmers.end() ){
            mean += binmer.second;
            }
        }
        return mean;
    }

    double _effectiveLength( TranscriptDataT* ts ) {
        double length = 0.0;
        for ( auto binmer : ts->binMers ) {
            length += _weight(binmer.first);
        }
        return length;
    }

    void _initialize( const std::vector<std::string> &transcriptFiles ) { 

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

        // Open up the transcript file for reading
        //const char *fnames[] = { transcriptFile.c_str() };
        // Create a jellyfish parser
        std::cerr << "here1\n";
        jellyfish::parse_read parser( fnames, fnames+numFnames, 1000);
        std::cerr << "here2\n";

        // So we can concisely identify each transcript
        TranscriptIDT transcriptIndex {0};

        // Holds the tasks to be processed
        boost::lockfree::fifo< TranscriptJob * > workQueue;
        // Holds the processed transcript data
        boost::lockfree::fifo< TranscriptDataT * > processedReadQueue;
        // Holds results that we will use to update the index map
        boost::lockfree::fifo< TranscriptKmerSet * > q;

        size_t numActors = 12;
        std::vector<std::thread> threads;

        _transcripts.resize( _transcriptGeneMap.numTranscripts(), nullptr );

        bool done {false};
        std::atomic<size_t> numRes {0};

        // Start the thread that will print the progress bar
        threads.push_back( std::thread( [&numRes, this] () {
            auto numJobs = this->_transcriptGeneMap.numTranscripts();
            boost::progress_display show_progress( numJobs );
            size_t lastCount = numRes;
            while ( lastCount < numJobs ) {
                auto diff = numRes - lastCount;
                if ( diff > 0 ) {
                    show_progress += diff;
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
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

                        // The set of kmers in this transcript
                        //auto ts = new TranscriptKmerSet {transcriptIndex, {}};

                        // We know how many kmers there will be, so allocate the space
                        size_t numKmers {readLen - this->_merLen + 1};
                        //std::get<1>(*ts).reserve(numKmers);
                        
                        // The binary representation of the kmers in this transcript
                        std::map<KmerIDT, KmerQuantityT> binMers;

                        size_t firstIndex{0};
                        // Iterate over the kmers
                        for ( auto offset : boost::irange( firstIndex, numKmers ) ) { // < numKmers ) {
                            auto mer = newSeq.substr( offset, this->_merLen );
                            // Only count and track kmers which should be considered
                            auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
                            if ( this->_considered(binMer) ) {                                
                                binMers[binMer] += 1.0;
                                //std::get<1>(*ts).push_back(binMer);
                                this->_transcriptsForKmer[binMer].push_back(transcriptIndex);
                            }
                        }
                        
                        auto procRead = new TranscriptInfo {
                            binMers, // the set of kmers
                            0.0, // initial desired mean
                            readLen // real transcript length
                        };
                        this->_transcripts[transcriptIndex] = procRead;
                    }
            }) );

        }

        // Wait for all of the threads to finish
        for ( auto& thread : threads ){ thread.join(); }

        numRes = 0;
        /*
        auto printProgressBar2 = [&numRes, this] () {
            auto numJobs = this->_transcriptsForKmer.size();
            boost::progress_display show_progress( numJobs );
            size_t lastCount = numRes;
            while ( lastCount < numJobs ) {
                auto diff = numRes - lastCount;
                if ( diff > 0 ) {
                    show_progress += diff;
                    lastCount += diff;
                }
                boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
            }
        };
        std::cerr << "Processing transcripts\n";
        tp.schedule(printProgressBar2);
        */

        tbb::task_scheduler_init tbb_init;

        std::cerr << "Determining gene-promiscuous kmers ... ";

        tbb::parallel_for_each( _transcriptsForKmer.begin(), _transcriptsForKmer.end(),
            [&numRes, this]( KmerMapT::value_type& kt ) {
                auto& transcripts = kt.second;
                std::sort(transcripts.begin(), transcripts.end());
                auto it = std::unique(transcripts.begin(), transcripts.end()); 
                transcripts.resize( it - transcripts.begin() );
                std::unordered_set<uint32_t> geneSet;
                for ( auto t : transcripts ){
                    geneSet.insert(this->_transcriptGeneMap.gene(t));
                }
                if ( geneSet.size() > 1 ){ this->_genePromiscuousKmers.insert(kt.first); }
                ++numRes;
            }
        );
        std::cerr << "done\n";

        std::cerr << "Computing initial coverage estimates ... ";

        tbb::parallel_for( size_t(0), size_t(_transcriptGeneMap.numTranscripts()),
        [this]( size_t tid ) {
            auto transcriptData = this->_transcripts[tid];
            for ( auto & kv : transcriptData->binMers ) {
                auto kmer = kv.first;
                if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                    auto count = kv.second;
                    kv.second = count * this->_readHash[kmer] * this->_weight(kmer);
                }
            }

            transcriptData->mean = this->_computeMean( transcriptData->binMers );
        }
        );


        std::cerr << "done\n";
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

        // For each kmer
        for ( auto & kmerTranscripts :  _transcriptsForKmer ) {
            auto kmer = kmerTranscripts.first;
            kmerList[idx] = kmer;
            ++idx;
        }
        
        tbb::task_scheduler_init tbb_init;

        for ( size_t iter = 0; iter < numIt; ++iter ) {

            auto reqNumJobs = _transcriptsForKmer.size();

            std::cerr << "iteraton: " << iter;

            globalError = 0.0;
            numJobs = 0;
            completedJobs = 0;

            auto pbthread = std::thread( 
                [&completedJobs, reqNumJobs]() -> bool {
                    auto prevNumJobs = 0;
                    boost::progress_display show_progress( reqNumJobs );
                    while ( prevNumJobs < reqNumJobs ) {
                        if ( prevNumJobs < completedJobs ) {
                            show_progress += completedJobs - prevNumJobs;
                        }
                        prevNumJobs = completedJobs.load();

                        boost::this_thread::sleep_for(boost::chrono::milliseconds(150));
                    }
                    return true;
                }
            );

            //tp.schedule( pbthread );


            tbb::parallel_for( size_t(0), size_t(kmerList.size()),
                
                [&kmerList, &completedJobs, this]( size_t kid ) {
                        // this thread claimed myJobID;
                        auto kmer = kmerList[kid];
                        if ( this->_genePromiscuousKmers.find(kmer) == this->_genePromiscuousKmers.end() ){
                        auto &transcripts = this->_transcriptsForKmer[kmer];
                        if ( transcripts.size() > 1 ) {

                            double totalMass = 0.0;
                            for ( auto tid : transcripts ) {
                                totalMass += this->_transcripts[tid]->mean;
                            }

                            double norm = 1.0 / totalMass;
                            for ( auto tid : transcripts ) {
                                this->_transcripts[tid]->binMers[kmer] =
                                    norm * this->_readHash[kmer] * this->_transcripts[tid]->mean;
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

            std::cerr << "\ncomputing new means ... ";
            // compute the new mean for each transcript
            tbb::parallel_for( size_t(0), size_t(_transcripts.size()),
                [this]( size_t tid ) -> void {
                        // this thread claimed myJobID;
                        auto ts = this->_transcripts[tid];
                        ts->mean = this->_computeMean( ts->binMers );
                }
            );
            std::cerr << "done\n";
        }


        std::ofstream ofile( outputFile );
        size_t index = 0;
        for ( auto ts : _transcripts ) {
            ofile << _transcriptGeneMap.transcriptName(index) << '\t' << ts->length << '\t' << _computeSum(ts->binMers) << "\n";
            ++index;
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