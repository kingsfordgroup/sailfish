#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "transcript_gene_map.hpp"
#include "genomic_feature.hpp"
#include "utils.hpp"

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


}