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


#ifndef __PERFECT_HASH_INDEX_HPP__
#define __PERFECT_HASH_INDEX_HPP__

#include <vector>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <memory>
#include <functional>

#include "boost/timer/timer.hpp"
#include "cmph.h"

//template <typename Deleter>
class PerfectHashIndex {
  typedef uint64_t Kmer;
  typedef uint32_t Count;
  typedef std::atomic<Kmer> AtomicKmer;
  typedef std::atomic<Count> AtomicCount;
  typedef std::function<void(cmph_t*)> Deleter;

  public:
   // We'll return this invalid id if a kmer is not found in our DB
   size_t INVALID = std::numeric_limits<size_t>::max();

   PerfectHashIndex( std::vector<Kmer>& kmers, std::unique_ptr<cmph_t, Deleter>& hash, uint32_t merSize ) : 
     kmers_(std::move(kmers)), hash_(std::move(hash)), merSize_( merSize ) {}

   PerfectHashIndex( PerfectHashIndex&& ph ) {
   	merSize_ = ph.merSize_;
   	hash_ = std::move(ph.hash_);
   	kmers_ = std::move(ph.kmers_);
   }

   void dumpToFile(const std::string& fname) {
   	FILE* out = fopen(fname.c_str(), "w");

   	// read the key set
    fwrite( reinterpret_cast<char*>(&merSize_), sizeof(merSize_), 1, out );
    size_t numCounts = kmers_.size();
    fwrite( reinterpret_cast<char*>(&numCounts), sizeof(size_t), 1, out );
    fwrite( reinterpret_cast<char*>(&kmers_[0]), sizeof(Kmer), numCounts, out );

    cmph_dump(hash_.get(), out); 
    fclose(out);

   }

   static PerfectHashIndex fromFile( const std::string& fname ) {
   	FILE* in = fopen(fname.c_str(),"r");

   	// read the key set
    uint32_t merSize;
    fread( reinterpret_cast<char*>(&merSize), sizeof(merSize), 1, in );
    size_t numCounts;
    fread( reinterpret_cast<char*>(&numCounts), sizeof(size_t), 1, in );
    std::vector<Kmer> kmers(numCounts, Kmer(0));
    fread( reinterpret_cast<char*>(&kmers[0]), sizeof(Kmer), numCounts, in );

    // read the hash
    std::unique_ptr<cmph_t, Deleter> hash( cmph_load(in), cmph_destroy );
    PerfectHashIndex index(kmers, hash, merSize);

    fclose(in);

    return index;
   }

   inline size_t index( uint64_t kmer ) {
   	char *key = (char*)(&kmer);
    unsigned int id = cmph_search(hash_.get(), key, sizeof(uint64_t));
    return (kmers_[id] == kmer) ? id : INVALID;
   }

   inline size_t numKeys() { return kmers_.size(); }

   bool verify() {
   	auto start = std::chrono::steady_clock::now();
   	for ( auto k : kmers_ ) { 
   		if( kmers_[index(k)] != k ) { return false; }
   	}
   	auto end = std::chrono::steady_clock::now();
   	auto ms = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
   	std::cerr << "verified: " << static_cast<double>(ms.count()) / kmers_.size() << " us / key\n";
   	return true;
   }

   inline uint32_t kmerLength() { return merSize_; }
   const std::vector<Kmer>& kmers() { return kmers_; }

   private:
   	std::vector<Kmer> kmers_;
   	std::unique_ptr<cmph_t, Deleter> hash_;
   	uint32_t merSize_;
};

#endif // __PERFECT_HASH_INDEX_HPP__