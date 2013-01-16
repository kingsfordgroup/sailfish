#ifndef COUNTDB_HPP
#define COUNTDB_HPP

#include <vector>
#include <atomic>
#include <fstream>
#include <limits>
#include <algorithm>

/**
*  This class provides low-overhead access to the counts of various
*  kmers in a hash-like format (though internally it is represented)
*  without hashing.
**/
class CountDB {
  typedef uint64_t Kmer;
  typedef uint32_t Count;
  typedef std::atomic<Kmer> AtomicKmer;
  typedef std::atomic<Count> AtomicCount;

  public:
   // We'll return this invalid id if a kmer is not found in our DB
   size_t INVALID = std::numeric_limits<size_t>::max();

   CountDB( std::vector<Kmer>& kmers, uint32_t merSize ) : 
           _kmers(kmers), _counts( std::vector< AtomicCount >( kmers.size() ) ), _merSize( merSize )
   {
    assert( std::is_sorted( kmers.begin(), kmers.end() ) );
   }

   /**
   *  This constructor populates this CountDB from file, reading the set of counts from fname.
   */
   CountDB( const std::string& fname, std::vector<Kmer>& kmers, uint32_t merSize ) :  
           _kmers(kmers), _merSize( merSize ) {

            assert( std::is_sorted( kmers.begin(), kmers.end() ) );
        std::ifstream counts(fname, std::ios::in | std::ios::binary );
        size_t numCounts{0};
        counts.read( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
        _counts = std::vector< AtomicCount >( numCounts );
        counts.read( reinterpret_cast<char*>(&_counts[0]), sizeof(_counts[0]) * numCounts );
        counts.close();

   }

   size_t id(Kmer k) {
    auto it = std::lower_bound( _kmers.begin(), _kmers.end(), k );
    if ( *it != k ) { return INVALID; }
    return std::distance(_kmers.begin(), it);
   }

   uint32_t operator[](uint64_t kmer) {
    auto idx = id(kmer);
    return (idx == INVALID) ? 0 : _counts[idx].load();
   }

   // increment the count for kmer 'k' by 'amt'
   // returns true if k existed in the database and false otherwise
   bool inc(Kmer k, uint32_t amt=1) {    
    auto it = std::lower_bound( _kmers.begin(), _kmers.end(), k );

    if ( *it != k ) { return false; }
    auto idx = std::distance(_kmers.begin(), it);
    
    _counts[idx] += amt;
    return true;
   }

   bool dumpToFile( const std::string& fname ){
    std::ofstream counts(fname, std::ios::out | std::ios::binary );
    size_t numCounts = _kmers.size();
    counts.write( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
    counts.write( reinterpret_cast<char*>(&_counts[0]), sizeof(_counts[0]) * numCounts );
    counts.close();
   }

  private:
    std::vector< Kmer >& _kmers;
    std::vector< AtomicCount > _counts;
    uint32_t _merSize;
};


#endif // COUNTDB_HPP