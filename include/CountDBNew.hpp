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


#ifndef COUNTDBNEW_HPP
#define COUNTDBNEW_HPP

#include <vector>
#include <atomic>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>

#include <sys/mman.h>

#include "tbb/concurrent_hash_map.h"
#include "PerfectHashIndex.hpp"

/**
*  This class provides low-overhead access to the counts of various
*  kmers in a hash-like format (though internally it may be represented)
*  without hashing.
**/
class CountDBNew {
  using Kmer = uint64_t;
  using Count = uint32_t;
  using AtomicCount = std::atomic<Count>;
  using Length = uint64_t;
  using AtomicLength = std::atomic<Length>;
  using AtomicLengthCount = std::atomic<Length>;

  public:
   // We'll return this invalid id if a kmer is not found in our DB
   size_t INVALID = std::numeric_limits<size_t>::max();

   CountDBNew( std::shared_ptr<PerfectHashIndex>& index ) :
      index_(index), counts_( std::vector< AtomicCount >( index->numKeys() ) ),
      length_(0), numLengths_(0) {}

   CountDBNew( CountDBNew&& other ) {
    counts_ = std::move(other.counts_);
    index_ = other.index_;
    length_ = other.length_.load();
    numLengths_ = other.numLengths_.load();
   }

   static CountDBNew fromFile( const std::string& fname, std::shared_ptr<PerfectHashIndex>& index ) {
    std::ifstream in(fname, std::ios::in | std::ios::binary );

    // Read in the total read length and # of reads
    uint64_t length = 0;
    uint64_t numLengths = 0;
    in.read(reinterpret_cast<char*>(&length), sizeof(length));
    in.read(reinterpret_cast<char*>(&numLengths), sizeof(numLengths));

    std::cerr << "read length = " << length << ", numLengths = " << numLengths << "\n";
    // Read in the count vector
    std::vector<AtomicCount> counts(index->numKeys());
    in.read( reinterpret_cast<char*>(&counts[0]), sizeof(AtomicCount) * index->numKeys() );
    in.close();

    CountDBNew cdb(index);
    cdb.counts_ = std::move(counts);
    cdb.length_ = length;
    cdb.numLengths_ = numLengths;
    return cdb;
   }

   inline void appendLength(uint32_t l, uint64_t count=1) {
    length_ += l;
    numLengths_ += count;
   }

   inline double averageLength() {
    return (numLengths_ > 0) ? (static_cast<double>(length_) / numLengths_) : 0.0;
   }

   inline Length totalLength() { return length_.load(); }
   inline Length numLengths() { return numLengths_.load(); }

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
   inline bool inc(Kmer k, uint32_t amt=1) {
    auto idx = index_->index(k);
    bool valid = (idx != INVALID);
    if (valid) { counts_[idx] += amt; }
    return valid;
   }

   inline void incAtIndex(std::vector<AtomicCount>::size_type idx, uint32_t amt=1) {
     counts_[idx] += amt;
   }

   void will_need(uint32_t threadIdx, uint32_t numThreads) {
     auto pageSize = sysconf(_SC_PAGESIZE);
     size_t numPages{0};

     auto entriesPerPage = pageSize / sizeof(AtomicCount);
     auto size = counts_.size();
     numPages = (sizeof(AtomicCount) * counts_.size()) / entriesPerPage;
     // number of pages that each thread should touch
     auto numPagesPerThread = numPages / numThreads;
     auto entriesPerThread = entriesPerPage * numPagesPerThread;
     // the page this thread starts touching
     auto start = entriesPerPage * threadIdx;
     for (size_t i = start; i < size; i += numThreads*entriesPerPage) {
      auto ci = counts_[i].load();
      counts_[i] = 0;
      counts_[i] = ci;
     }
   }

   bool dumpCountsToFile( const std::string& fname ) {
    std::ofstream counts(fname, std::ios::out | std::ios::binary );
    uint64_t length = length_.load();
    uint64_t numLengths = numLengths_.load();
    counts.write(reinterpret_cast<char*>(&length), sizeof(length));
    counts.write(reinterpret_cast<char*>(&numLengths), sizeof(numLengths));
    size_t numCounts = counts_.size();
    counts.write( reinterpret_cast<char*>(&counts_[0]), sizeof(counts_[0]) * numCounts );
    counts.close();
    return true;
   }

   inline uint32_t kmerLength() { return index_->kmerLength(); }
   const std::vector<Kmer>& kmers() { return index_->kmers(); }
  private:
    std::shared_ptr<PerfectHashIndex> index_;
    std::vector< AtomicCount > counts_;
    AtomicLength length_;
    AtomicLengthCount numLengths_;
};


#endif // COUNTDBNEW_HPP
