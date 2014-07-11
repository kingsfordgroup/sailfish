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
#include <boost/filesystem.hpp>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>

#include "gff.h"

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "jellyfish/mer_dna.hpp"
#include "TranscriptGeneMap.hpp"
#include "GenomicFeature.hpp"
#include "SailfishUtils.hpp"

namespace sailfish {
namespace utils {
using std::string;
using NameVector = std::vector<string>;
using IndexVector = std::vector<size_t>;
using KmerVector = std::vector<uint64_t>;

/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatStringNew(std::string& fmt) {
	using std::vector;
	using std::string;
	using std::map;
	using std::stringstream;

    map<string, LibraryFormat> formatMap = {
        {"IU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U)},
        {"ISF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA)},
        {"ISR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS)},
        {"OU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::U)},
        {"OSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA)},
        {"OSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS)},
        {"MU", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::U)},
        {"MSF", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S)},
        {"MSR", LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A)},
        {"U", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U)},
        {"SF", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S)},
        {"SR", LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::A)}};

	// inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
	// first convert the string to upper-case
	for (auto& c : fmt) { c = std::toupper(c); }


    auto libFmtIt = formatMap.find(fmt);

	if (libFmtIt == formatMap.end()) {
		stringstream errstr;
		errstr << "unknown library format string : " << fmt;
		throw std::invalid_argument(errstr.str());
	}

    return libFmtIt->second;
}

/**
 * Parses a set of __ordered__ command line options and extracts the relevant
 * read libraries from them.
 */
std::vector<ReadLibrary> extractReadLibraries(boost::program_options::parsed_options& orderedOptions) {
	// The current (default) format for paired end data
	LibraryFormat peFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::U);
	// The current (default) format for single end data
	LibraryFormat seFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::U);

	std::vector<ReadLibrary> peLibs{ReadLibrary(peFormat)};
	std::vector<ReadLibrary> seLibs{ReadLibrary(seFormat)};
	for (auto& opt : orderedOptions.options) {
		// Update the library type
		if (opt.string_key == "libtype") {
			auto libFmt = parseLibraryFormatStringNew(opt.value[0]);
			if (libFmt.type == ReadType::PAIRED_END) {
				peFormat = libFmt;
				peLibs.emplace_back(libFmt);
			} else {
				seFormat = libFmt;
				seLibs.emplace_back(libFmt);
			}
		}
		if (opt.string_key == "mates1") {
			std::cerr << "mates1\n";
			peLibs.back().addMates1(opt.value);
		}
		if (opt.string_key == "mates2") {
			std::cerr << "mates2\n";
			peLibs.back().addMates2(opt.value);
		}
		if (opt.string_key == "unmated_reads") {
			std::cerr << "unmated\n";
			seLibs.back().addUnmated(opt.value);
		}
	}

	std::vector<ReadLibrary> libs;
	libs.reserve(peLibs.size() + seLibs.size());
	for (auto& lib : boost::range::join(seLibs, peLibs)) {
		if (lib.format().type == ReadType::SINGLE_END) {
			if (lib.unmated().size() == 0) {
				std::cerr << "skipping single-end library w/ no reads\n";
				continue;
			}
		} else if (lib.format().type == ReadType::PAIRED_END) {
			if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
				std::cerr << "skipping paired-end library w/ no reads\n";
				continue;
			}
		}
		libs.push_back(lib);
	}
	std::cerr << "there are " << libs.size() << " libs\n";
	return libs;
}



/**
 * This function parses the library format string that specifies the format in which
 * the reads are to be expected.
 */
LibraryFormat parseLibraryFormatString(std::string& fmt) {
    using std::vector;
    using std::string;
    using std::map;
    using std::stringstream;

    // inspired by http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c

    // first convert the string to upper-case
    for (auto& c : fmt) { c = std::toupper(c); }
    // split on the delimiter ':', and put the key, value (k=v) pairs into a map
    stringstream ss(fmt);
    string item;
    map<string, string> kvmap;
    while (std::getline(ss, item, ':')) {
        auto splitPos = item.find('=', 0);
        string key{item.substr(0, splitPos)};
        string value{item.substr(splitPos+1)};
        kvmap[key] = value;
    }

    map<string, ReadType> readType = {{"SE", ReadType::SINGLE_END}, {"PE", ReadType::PAIRED_END}};
    map<string, ReadOrientation> orientationType = {{">>", ReadOrientation::SAME},
                                           {"<>", ReadOrientation::AWAY},
                                           {"><", ReadOrientation::TOWARD},
                                           {"*", ReadOrientation::NONE}};
    map<string, ReadStrandedness> strandType = {{"SA", ReadStrandedness::SA},
                                    {"AS", ReadStrandedness::AS},
                                    {"A", ReadStrandedness::A},
                                    {"S", ReadStrandedness::S},
                                    {"U", ReadStrandedness::U}};
    auto it = kvmap.find("T");
    string typeStr = "";
    if (it != kvmap.end()) {
        typeStr = it->second;
    } else {
        it = kvmap.find("TYPE");
        if (it != kvmap.end()) {
            typeStr = it->second;
        }
    }

    if (typeStr != "SE" and typeStr != "PE") {
        string e = typeStr + " is not a valid read type; must be one of {SE, PE}";
        throw std::invalid_argument(e);
    }

    ReadType type = (typeStr == "SE") ? ReadType::SINGLE_END : ReadType::PAIRED_END;
    ReadOrientation orientation = (type == ReadType::SINGLE_END) ? ReadOrientation::NONE : ReadOrientation::TOWARD;
    ReadStrandedness strandedness{ReadStrandedness::U};
    // Construct the LibraryFormat class from the key, value map
    for (auto& kv : kvmap) {
        auto& k = kv.first; auto& v = kv.second;
        if (k == "O" or k == "ORIENTATION") {
            auto it = orientationType.find(v);
            if (it != orientationType.end()) { orientation = orientationType[it->first]; } else {
                string e = v + " is not a valid orientation type; must be one of {>>, <>, ><}";
                throw std::invalid_argument(e);
            }

        }
        if (k == "S" or k == "STRAND") {
            auto it = strandType.find(v);
            if (it != strandType.end()) { strandedness = strandType[it->first]; } else {
                string e = v + " is not a valid strand type; must be one of {SA, AS, S, A, U}";
                throw std::invalid_argument(e);
            }
        }

    }
    LibraryFormat lf(type, orientation, strandedness);
    return lf;
}



uint64_t encode(uint64_t tid, uint64_t offset) {
    uint64_t res = (((tid & 0xFFFFFFFF) << 32) | (offset & 0xFFFFFFFF));
    return res;
}

uint32_t transcript(uint64_t enc) {
    uint32_t t = (enc & 0xFFFFFFFF00000000) >> 32;
    return t;
}

uint32_t offset(uint64_t enc) {
    uint32_t o = enc & 0xFFFFFFFF;
    return o;
}

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


TranscriptGeneMap transcriptGeneMapFromGTF(const std::string& fname, std::string key) {

    using std::unordered_set;
    using std::unordered_map;
    using std::vector;
    using std::tuple;
    using std::string;
    using std::get;

    // Use GffReader to read the file
    GffReader reader(const_cast<char*>(fname.c_str()));
    // Remember the optional attributes
    reader.readAll(true);

    struct TranscriptKeyPair {
        const char* transcript_id;
        const char* key;
        TranscriptKeyPair(const char* t, const char* k) :
            transcript_id(t), key(k) {}
    };

    // The user can group transcripts by gene_id, gene_name, or
    // an optinal attribute that they provide as a string.
    enum class TranscriptKey { GENE_ID, GENE_NAME, DYNAMIC };

    // Select the proper attribute by which to group
    TranscriptKey tkey = TranscriptKey::GENE_ID;

    if (key == "gene_id") {
    } else if (key == "gene_name") {
        tkey = TranscriptKey::GENE_NAME;
    } else {
        tkey = TranscriptKey::DYNAMIC;
    }

    // Iterate over all transcript features and build the
    // transcript <-> key vector.
    auto nfeat = reader.gflst.Count();
    std::vector<TranscriptKeyPair> feats;
    for (int i=0; i < nfeat; ++i) {
        auto f = reader.gflst[i];
        if (f->isTranscript()) {
            const char* keyStr;
            switch (tkey) {
                case TranscriptKey::GENE_ID:
                    keyStr = f->getGeneID();
                    break;
                case TranscriptKey::GENE_NAME:
                    keyStr = f->getGeneName();
                    break;
                case TranscriptKey::DYNAMIC:
                    keyStr = f->getAttr(key.c_str());
                    break;
            }
            feats.emplace_back(f->getID(), keyStr);
        }
    }

    // Given the transcript <-> key vector, build the
    // TranscriptGeneMap.

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
    []( const TranscriptKeyPair & a, const TranscriptKeyPair & b) -> bool {
        return std::strcmp(a.transcript_id, b.transcript_id) < 0;
    } );

    std::string currentTranscript = "";
    for ( auto & feat : feats ) {

        std::string gene(feat.key);
        std::string transcript(feat.transcript_id);

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


TranscriptGeneMap readTranscriptToGeneMap( std::ifstream &ifile ) {

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
    using std::vector;
    using stream_manager = jellyfish::stream_manager<char**>;
    using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;
    namespace bfs = boost::filesystem;

    NameVector transcriptNames;
    NameVector geneNames {"gene"};

    vector<bfs::path> paths{transcriptsFile};

    // Create a jellyfish parser
    const int concurrentFile{1};
    char** fnames = new char*[1];
    fnames[0] = const_cast<char*>(transcriptsFile.c_str());
    stream_manager streams(fnames, fnames + 1, concurrentFile);

    size_t maxReadGroupSize{100};
    sequence_parser parser(4, maxReadGroupSize, concurrentFile, streams);

    // while there are transcripts left to process
    while (true) {
        sequence_parser::job j(parser);
        // If this job is empty, then we're done
        if (j.is_empty()) { break; }

        for (size_t i=0; i < j->nb_filled; ++i) {
            // The transcript name
            std::string fullHeader(j->data[i].header);
            std::string header = fullHeader.substr(0, fullHeader.find(' '));
            transcriptNames.emplace_back(header);
        }
    }

    // Sort the transcript names
    std::sort(transcriptNames.begin(), transcriptNames.end());

    // Since we have no real gene groupings, the t2g vector is trivial,
    // everything maps to gene 0.
    IndexVector t2g(transcriptNames.size(), 0);

    return TranscriptGeneMap(transcriptNames, geneNames, t2g);
}

class ExpressionRecord {
    public:
        ExpressionRecord() : name(""), length(0), quants(std::vector<double>()) {}
        std::string name;
        unsigned long length;
        std::vector<double> quants;
};

ExpressionRecord parseExpressionRecord(std::string& l) {
    // name
    // length
    // estimates
    ExpressionRecord rec;
    boost::char_separator<char> sep(" \t");
    boost::tokenizer<boost::char_separator<char>> tokens(l, sep);

    size_t idx{0};
    for (auto& tok : tokens) {
        switch (idx) {
            case 0:
                rec.name = tok;
                break;
            case 1:
                rec.length = stoul(tok);
                break;
            default:
                rec.quants.push_back(stod(tok));
        }
        ++idx;
    }

    return rec;
}

void aggregateEstimatesToGeneLevel(TranscriptGeneMap& tgm, boost::filesystem::path& inputPath) {

    std::ifstream inputFile(inputPath.string());
    std::string l;

    if (!inputFile.is_open()) {
        perror("Error reading file");
    }

    std::vector<std::string> comments;
    std::vector<ExpressionRecord> geneExpressions;
    geneExpressions.resize(tgm.numGenes());
    size_t numExpressionColumns{0};
    while (std::getline(inputFile, l)) {
        // If this is a comment line
        if (l.front() == '#') {
            comments.emplace_back(l);
        } else {
            auto expRec = parseExpressionRecord(l);
            numExpressionColumns = expRec.quants.size();
            auto tid = tgm.findTranscriptID(expRec.name);
            auto& geneExpressionRecord = geneExpressions[tgm.gene(tid)];
            size_t currLen = geneExpressionRecord.length;
            geneExpressionRecord.length = std::max(currLen, expRec.length);

            if (numExpressionColumns > geneExpressionRecord.quants.size()) {
                geneExpressionRecord.quants.resize(numExpressionColumns);
            }

            for (size_t ce = 0; ce < numExpressionColumns; ++ce) {
                geneExpressionRecord.quants[ce] += expRec.quants[ce];
            }
        }
    }
    inputFile.close();


    boost::filesystem::path outputFilePath(inputPath);
    outputFilePath.replace_extension(".genes.sf");
    std::ofstream outputFile(outputFilePath.string());
    for (auto& c : comments) {
        outputFile << c << "\n";
    }

    for (size_t gid = 0; gid < tgm.numGenes(); ++gid) {
        std::string geneName = tgm.nameFromGeneID(gid);
        auto& expRec = geneExpressions[gid];
        outputFile << geneName << '\t' << expRec.length;
        for (size_t ec = 0; ec < numExpressionColumns; ++ec) {
            outputFile << '\t' << expRec.quants[ec];
        }
        outputFile << '\n';
    }
    outputFile.close();
}

}
}
