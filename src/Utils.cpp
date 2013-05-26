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


#include <boost/thread/thread.hpp>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/parse_dna.hpp>

#include "TranscriptGeneMap.hpp"
#include "GenomicFeature.hpp"
#include "Utils.hpp"

namespace utils {
using std::string;
typedef std::vector<string> NameVector;
typedef std::vector<size_t> IndexVector;
typedef std::vector<uint64_t> KmerVector;

size_t numberOfReadsInFastaFile(const std::string& fname) {
    constexpr size_t bufferSize = 16184;
    char buffer[bufferSize];
    std::ifstream ifile(fname, std::ifstream::in);
    ifile.rdbuf()->pubsetbuf(buffer, bufferSize);

    size_t numReads = 0;
    std::string s;
    while (ifile >> s) { if (s.front() == '>') { ++numReads; } }

    ifile.close();

    return numReads;
}

bool readKmerOrder( const std::string& fname, std::vector<uint64_t>& kmers ) {

  std::ifstream mlist(fname, std::ios::in | std::ios::binary);
  // Get the number of kmers from file
  size_t numKmers{0};
  mlist.read( reinterpret_cast<char*>( &numKmers ), sizeof( size_t ) );
  
  // Resize the array that will hold the sorted kmers
  kmers.resize(numKmers, 0);
  mlist.read( reinterpret_cast<char*>( &kmers[0] ), sizeof( uint64_t) * kmers.size() ); 

  mlist.close();

  return true;
}

template <template<typename> class S, typename T>
bool overlap( const S<T> &a, const S<T> &b ) {
    // Query from the smaller set to the larger set
    if ( a.size() <= b.size() ) {
        for ( auto & ae : a ) {
            if (b.find(ae) != b.end()) {
                return true;
            }
        }
    } else {
        for ( auto & be : b ) {
            if (a.find(be) != b.end()) {
                return true;
            }
        }
    }
    // If nothing from the smaller set is in the larger set, then they don't overlap
    return false;
}

TranscriptGeneMap readTranscriptToGeneMap( std::ifstream &ifile ) {

    using std::unordered_set;
    using std::unordered_map;
    using std::vector;
    using std::tuple;
    using std::string;
    using std::get;

    typedef tuple<string, size_t> NameID;

    IndexVector t2g;
    NameVector transcriptNames;
    NameVector geneNames;

    // holds the transcript name ID mapping
    vector<NameID> transcripts;
    // holds the mapping from transcript ID to gene ID
    IndexVector t2gUnordered;
    // holds the set of gene IDs
    unordered_map<string, size_t> geneNameToID;

    // To read the input and assign ids
    size_t transcriptCounter = 0;
    size_t geneCounter = 0;
    string transcript;
    string gene;

    while ( ifile >> transcript >> gene ) {
        // The transcript and it's ID
        transcripts.push_back( make_tuple(transcript, transcriptCounter) );

        auto geneIt = geneNameToID.find(gene);
        size_t geneID = 0;

        if ( geneIt == geneNameToID.end() ) {
            // If we haven't seen this gene yet, give it a new ID
            geneNameToID[gene] = geneCounter;
            geneID = geneCounter;
            geneNames.push_back(gene);
            ++geneCounter;
        } else {
            // Otherwise lookup the ID
            geneID = geneIt->second;
        }

        // Map the transcript to the gene in terms of their IDs
        t2gUnordered.push_back(geneID);

        ++transcriptCounter;
    }

    std::sort( transcripts.begin(), transcripts.end(),
               []( const NameID & a, const NameID & b) -> bool { return get<0>(a) < get<0>(b); } );

    // Resize these vectors for fast access
    transcriptNames.resize(t2gUnordered.size());
    t2g.resize(t2gUnordered.size());

    for ( size_t newID = 0; newID < transcripts.size(); ++newID ) {
        // For each transcript, map it to the appropriate gene
        string oldName; size_t oldID;
        std::tie(oldName, oldID) = transcripts[newID];
        t2g[newID] = t2gUnordered[oldID];
        transcriptNames[newID] = oldName;
    }

    return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}


TranscriptGeneMap transcriptToGeneMapFromFasta( const std::string& transcriptsFile ) {
    typedef std::vector<string> NameVector;
    typedef std::vector<size_t> IndexVector;

    NameVector transcriptNames;
    NameVector geneNames {"gene"};

    char* fnames[1]  = { const_cast<char*>(transcriptsFile.c_str()) };
    // Create a jellyfish parser
    jellyfish::parse_read parser( fnames, fnames+1, 1000);

    // Each thread gets it's own stream
    jellyfish::parse_read::thread stream = parser.new_thread();
    jellyfish::parse_read::read_t* read;

    // while there are transcripts left to process
    while ( (read = stream.next_read()) ) { 
      // The transcript name
      std::string fullHeader(read->header, read->hlen);
      std::string header = fullHeader.substr(0, fullHeader.find(' '));
      transcriptNames.emplace_back(header);
    }

    // Sort the transcript names
    std::sort(transcriptNames.begin(), transcriptNames.end());
    
    // Since we have no real gene groupings, the t2g vector is trivial,
    // everything maps to gene 0.
    IndexVector t2g(transcriptNames.size(), 0);

    return TranscriptGeneMap(transcriptNames, geneNames, t2g);

}

}