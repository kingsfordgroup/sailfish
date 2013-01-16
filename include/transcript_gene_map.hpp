#ifndef TRANSCRIPT_GENE_MAP_HPP
#define TRANSCRIPT_GENE_MAP_HPP

#include <algorithm>
#include <vector>

class TranscriptGeneMap {
    typedef size_t Index;
    typedef size_t Size;
    typedef std::vector<std::string> NameVector;
    typedef std::vector<size_t> IndexVector;
    typedef std::vector<std::vector<size_t>> IndexVectorList;

private:
    NameVector _transcriptNames;
    NameVector _geneNames;
    IndexVector _transcriptsToGenes;
    IndexVectorList _genesToTranscripts;
    bool _haveReverseMap;

    void _computeReverseMap() {

        _genesToTranscripts.resize( _geneNames.size(), {});

        size_t geneID;
        size_t transcriptID = 0;
        size_t maxNumTrans = 0;
        std::string maxGene;
        for ( size_t transcriptID = 0; transcriptID < _transcriptsToGenes.size(); ++transcriptID ) {
            _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].push_back( transcriptID );
            if ( maxNumTrans < _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].size() ) {
                maxNumTrans = _genesToTranscripts[ _transcriptsToGenes[transcriptID] ].size();
                maxGene = _transcriptsToGenes[transcriptID];
            }
        }
        std::cerr << "max # of transcripts in a gene was " << maxNumTrans << " in gene " << maxGene << "\n";
    }

public:
    TranscriptGeneMap( const NameVector &transcriptNames,
                       const NameVector &geneNames,
                       const IndexVector &transcriptsToGenes ) :
        _transcriptNames(transcriptNames), _geneNames(geneNames),
        _transcriptsToGenes(transcriptsToGenes), _haveReverseMap(false) {}

    Index INVALID { std::numeric_limits<Index>::max() };

    Index findTranscriptID( const std::string &tname ) {
        using std::distance;
        using std::lower_bound;
        auto it = lower_bound( _transcriptNames.begin(), _transcriptNames.end(), tname );
        return ( it == _transcriptNames.end() ) ? INVALID : ( distance(_transcriptNames.begin(), it) );
    }

    Size numTranscripts() {
        return _transcriptNames.size();
    }
    Size numGenes() {
        return _geneNames.size();
    }

    bool needReverse() {
        if ( _haveReverseMap ) {
            return false;
        } else {
            _computeReverseMap();
            return true;
        }
    }

    const IndexVector &transcriptsForGene( Index geneID ) {
        return _genesToTranscripts[geneID];
    }

    inline std::string nameFromGeneID( Index geneID ) {
        return _geneNames[geneID];
    }
    inline Index gene( Index transcriptID ) {
        return _transcriptsToGenes[transcriptID];
    }
    inline std::string geneName( Index transcriptID ) {
        return _geneNames[_transcriptsToGenes[transcriptID]];
    }
    inline std::string transcriptName( Index transcriptID ) {
        return _transcriptNames[transcriptID];
    }
};

#endif // TRANSCRIPT_GENE_MAP_HPP