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


#ifndef SAILFISH_UTILS_HPP
#define SAILFISH_UTILS_HPP

#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "TranscriptGeneMap.hpp"
#include "GenomicFeature.hpp"

namespace sailfish{
namespace utils {
using std::string;
using NameVector = std::vector<string>;
using IndexVector = std::vector<size_t>;
using KmerVector = std::vector<uint64_t>;

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

    using NameID = tuple<string, size_t>;

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

TranscriptGeneMap transcriptToGeneMapFromFasta( const std::string& transcriptsFile );

}
}
#endif // UTILS_HPP
