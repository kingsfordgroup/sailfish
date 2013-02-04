#ifndef COUNTDBNEW_HPP
#define COUNTDBNEW_HPP

#include <vector>
#include <atomic>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>

#include "PerfectHashIndex.hpp"

/**
*  This class provides low-overhead access to the counts of various
*  kmers in a hash-like format (though internally it is represented)
*  without hashing.
**/
class CountDBNew {
  typedef uint64_t Kmer;
  typedef uint32_t Count;
  typedef std::atomic<Count> AtomicCount;

  public:
   // We'll return this invalid id if a kmer is not found in our DB
   size_t INVALID = std::numeric_limits<size_t>::max();

   CountDBNew( std::shared_ptr<PerfectHashIndex>& index ) : 
      index_(index), counts_( std::vector< AtomicCount >( index->numKeys() ) ) {}

   CountDBNew( CountDBNew&& other ) {
    counts_ = std::move(other.counts_);
    index_ = other.index_;
   }

   static CountDBNew fromFile( const std::string& fname, std::shared_ptr<PerfectHashIndex>& index ) {
    std::ifstream in(fname, std::ios::in | std::ios::binary );
    std::vector<AtomicCount> counts(index->numKeys());
    in.read( reinterpret_cast<char*>(&counts[0]), sizeof(AtomicCount) * index->numKeys() );
    in.close();
    CountDBNew cdb(index);
    cdb.counts_ = std::move(counts);
    return cdb;
   }


   inline size_t id(Kmer k) { return index_->index(k); }

   uint32_t operator[](uint64_t kmer) {
    auto idx = id(kmer);
    return (idx == INVALID) ? 0 : counts_[idx].load();
   }

   uint32_t atIndex(size_t idx) {
      return (idx == INVALID) ? 0 : counts_[idx].load();
   }

   std::vector<AtomicCount>::size_type size() { return counts_.size(); }

   // increment the count for kmer 'k' by 'amt'
   // returns true if k existed in the database and false otherwise
   bool inc(Kmer k, uint32_t amt=1) {
    auto idx = index_->index(k);
    bool valid = (idx != INVALID);
    if (valid) { counts_[idx] += amt; }
    return valid;
   }

   bool dumpCountsToFile( const std::string& fname ) {
    std::ofstream counts(fname, std::ios::out | std::ios::binary );
    size_t numCounts = counts_.size();
    //counts.write( reinterpret_cast<char*>(&_merSize), sizeof(_merSize) );
    //counts.write( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
    counts.write( reinterpret_cast<char*>(&counts_[0]), sizeof(counts_[0]) * numCounts );
    counts.close();
   }

   inline uint32_t kmerLength() { return index_->kmerLength(); }
   const std::vector<Kmer>& kmers() { return index_->kmers(); }
  private:
    std::shared_ptr<PerfectHashIndex> index_;
    std::vector< AtomicCount > counts_;
};


#endif // COUNTDBNEW_HPP