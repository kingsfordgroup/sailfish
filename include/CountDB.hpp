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


#ifndef COUNTDB_HPP
#define COUNTDB_HPP

#include <vector>
#include <atomic>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>

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

   CountDB( std::shared_ptr<std::vector<Kmer>>& kmers, uint32_t merSize ) : 
           _kmers(kmers), _counts( std::vector< AtomicCount >( kmers->size() ) ), _merSize( merSize )
   {
    assert( std::is_sorted( kmers->begin(), kmers->end() ) );
   }

   static CountDB fromFile( const std::string& fname ) {
    std::ifstream in(fname, std::ios::in | std::ios::binary );
    uint32_t merSize;
    in.read( reinterpret_cast<char*>(&merSize), sizeof(merSize) );
    size_t numCounts;
    in.read( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
    auto kmers = std::make_shared<std::vector<Kmer>>(numCounts, Kmer(0));
    std::vector<AtomicCount> counts(numCounts);
    in.read( reinterpret_cast<char*>(&((*kmers)[0])), sizeof(Kmer) * numCounts );
    in.read( reinterpret_cast<char*>(&counts[0]), sizeof(AtomicCount) * numCounts );
    in.close();

    CountDB index(kmers, merSize);
    index._counts = std::move(counts);
    return index;
   }

   /**
   *  This constructor populates this CountDB from file, reading the set of counts from fname.
   */
   CountDB( const std::string& fname, std::shared_ptr<std::vector<Kmer>>& kmers, uint32_t merSize ) :  
     _kmers(kmers), _merSize( merSize ) {

      std::cerr << "checking that kmers are sorted . . . ";
      assert( std::is_sorted( kmers->begin(), kmers->end() ) );
      std::cerr << "done\n";


      std::ifstream counts(fname, std::ios::in | std::ios::binary );
      uint32_t fileMerSize{0};
      counts.read( reinterpret_cast<char*>(&fileMerSize), sizeof(fileMerSize) );
      std::cerr << "checking that the mersizes are the same . . .";
      assert( fileMerSize == merSize );
      std::cerr << "done\n";

      size_t numCounts{0};
      
      counts.read( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
      
      std::cerr << "checking that the counts are the same . . .";
      assert( numCounts == kmers->size() );
      std::cerr << "done\n";

      _counts = std::vector< AtomicCount >( numCounts );
      std::cerr << "reading counts from file . . .";
      counts.read( reinterpret_cast<char*>(&_counts[0]), sizeof(_counts[0]) * numCounts );
      std::cerr << "done\n";
      counts.close();
   }

   size_t id(Kmer k) {
    auto it = std::lower_bound( _kmers->begin(), _kmers->end(), k );
    if ( *it != k ) { return INVALID; }
    return std::distance(_kmers->begin(), it);
   }

   uint32_t operator[](uint64_t kmer) {
    auto idx = id(kmer);
    return (idx == INVALID) ? 0 : _counts[idx].load();
   }

   uint32_t atIndex(size_t idx) {
      return (idx == INVALID) ? 0 : _counts[idx].load();
   }

   std::vector<AtomicCount>::size_type size() { return _counts.size(); }

   // increment the count for kmer 'k' by 'amt'
   // returns true if k existed in the database and false otherwise
   bool inc(Kmer k, uint32_t amt=1) {    
    auto it = std::lower_bound( _kmers->begin(), _kmers->end(), k );

    if ( *it != k ) { return false; }
    auto idx = std::distance(_kmers->begin(), it);
    
    _counts[idx] += amt;
    return true;
   }

   bool dumpCountsToFile( const std::string& fname ) {
    std::ofstream counts(fname, std::ios::out | std::ios::binary );
    size_t numCounts = _kmers->size();
    counts.write( reinterpret_cast<char*>(&_merSize), sizeof(_merSize) );
    counts.write( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
    counts.write( reinterpret_cast<char*>(&_counts[0]), sizeof(_counts[0]) * numCounts );
    counts.close();
   }

   bool dumpToFile( const std::string& fname ){
    std::ofstream out(fname, std::ios::out | std::ios::binary );
    size_t numCounts = _kmers->size();
    out.write( reinterpret_cast<char*>(&_merSize), sizeof(_merSize) );
    out.write( reinterpret_cast<char*>(&numCounts), sizeof(size_t) );
    out.write( reinterpret_cast<char*>(&((*_kmers)[0])), sizeof((*_kmers)[0]) * numCounts );
    out.write( reinterpret_cast<char*>(&_counts[0]), sizeof(_counts[0]) * numCounts );
    out.close();
   }

   std::shared_ptr<std::vector< Kmer >>& indexKmers() { return _kmers; }
   uint32_t kmerLength() { return _merSize; }

  private:
    std::shared_ptr<std::vector< Kmer >> _kmers;
    std::vector< AtomicCount > _counts;
    uint32_t _merSize;
};


#endif // COUNTDB_HPP