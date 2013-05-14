#ifndef SAILFISH_UTILS_HPP
#define SAILFISH_UTILS_HPP

#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "transcript_gene_map.hpp"
#include "genomic_feature.hpp"

namespace utils {
using std::string;
typedef std::vector<string> NameVector;
typedef std::vector<size_t> IndexVector;
typedef std::vector<uint64_t> KmerVector;

size_t numberOfReadsInFastaFile(const std::string& fname);

bool readKmerOrder( const std::string& fname, std::vector<uint64_t>& kmers );

template <template<typename> class S, typename T>
bool overlap( const S<T> &a, const S<T> &b );

template< typename T >
TranscriptGeneMap transcriptToGeneMapFromFeatures( std::vector<GenomicFeature<T>> &feats ) {
    using std::unordered_set;
    using std::unordered_map;
    using std::vector;
    using std::tuple;
    using std::string;
    using std::get;

    typedef tuple<string, size_t> NameID;
    typedef std::vector<string> NameVector;
    typedef std::vector<size_t> IndexVector;

    IndexVector t2g;
    NameVector transcriptNames;
    NameVector geneNames;

    // holds the mapping from transcript ID to gene ID
    IndexVector t2gUnordered;
    // holds the set of gene IDs
    unordered_map<string, size_t> geneNameToID;

    // To read the input and assign ids
    size_t transcriptCounter = 0;
    size_t geneCounter = 0;
    string transcript;
    string gene;

    std::sort( feats.begin(), feats.end(),
    []( const GenomicFeature<T> & a, const GenomicFeature<T> & b) -> bool {
        return a.sattr.transcript_id < b.sattr.transcript_id;
    } );

    std::string currentTranscript = "";
    for ( auto & feat : feats ) {

        auto &gene = feat.sattr.gene_id;
        auto &transcript = feat.sattr.transcript_id;

        if ( transcript != currentTranscript ) {
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

            transcriptNames.push_back(transcript);
            t2g.push_back(geneID);

            //++transcriptID;
            currentTranscript = transcript;
        }

    }

    return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}

TranscriptGeneMap readTranscriptToGeneMap( std::ifstream &ifile );
}

#endif // UTILS_HPP
