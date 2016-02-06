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
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/join.hpp>

#include "spdlog/spdlog.h"

#include "gff.h"

#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "TranscriptGeneMap.hpp"
#include "SailfishUtils.hpp"
#include "ReadExperiment.hpp"

//S_AYUSH_CODE
#include "ReadKmerDist.hpp"
#include "UtilityFunctions.hpp"
//T_AYUSH_CODE

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
        LibraryFormat parseLibraryFormatString(std::string& fmt) {
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
                if (opt.string_key == "libType") {
                    auto libFmt = parseLibraryFormatString(opt.value[0]);
                    if (libFmt.type == ReadType::PAIRED_END) {
                        peFormat = libFmt;
                        peLibs.emplace_back(libFmt);
                    } else {
                        seFormat = libFmt;
                        seLibs.emplace_back(libFmt);
                    }
                }
                if (opt.string_key == "mates1") {
                    peLibs.back().addMates1(opt.value);
                }
                if (opt.string_key == "mates2") {
                    peLibs.back().addMates2(opt.value);
                }
                if (opt.string_key == "unmatedReads") {
                    seLibs.back().addUnmated(opt.value);
                }
            }

            std::vector<ReadLibrary> libs;
            libs.reserve(peLibs.size() + seLibs.size());
            for (auto& lib : boost::range::join(seLibs, peLibs)) {
                if (lib.format().type == ReadType::SINGLE_END) {
                    if (lib.unmated().size() == 0) {
                        // Didn't use default single end library type
                        continue;
                    }
                } else if (lib.format().type == ReadType::PAIRED_END) {
                    if (lib.mates1().size() == 0 or lib.mates2().size() == 0) {
                        // Didn't use default paired-end library type
                        continue;
                    }
                }
                libs.push_back(lib);
            }
            size_t numLibs = libs.size();
            std::cerr << "there " << ((numLibs > 1) ? "are " : "is ") << libs.size() << ((numLibs > 1) ? " libs\n" : " lib\n");
            return libs;
        }


        // for single end reads or orphans
        bool compatibleHit(LibraryFormat expected,
                           int32_t start, bool isForward, MateStatus ms) {
            auto expectedStrand = expected.strandedness;
            switch (ms) {
                case MateStatus::SINGLE_END:
                    if (isForward) { // U, SF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    } else { // U, SR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    }
                    break;
                case MateStatus::PAIRED_END_LEFT:
                    // "M"atching or same orientation is a special case
                    if (expected.orientation == ReadOrientation::SAME) {
                        return (expectedStrand == ReadStrandedness::U
                                or
                                (expectedStrand == ReadStrandedness::S and isForward)
                                or
                                (expectedStrand == ReadStrandedness::A and !isForward));
                    } else if (isForward) { // IU, ISF, OU, OSF, MU, MSF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    } else { // IU, ISR, OU, OSR, MU, MSR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    }
                    break;
                case MateStatus::PAIRED_END_RIGHT:
                    // "M"atching or same orientation is a special case
                    if (expected.orientation == ReadOrientation::SAME) {
                        return (expectedStrand == ReadStrandedness::U
                                or
                                (expectedStrand == ReadStrandedness::S and isForward)
                                or
                                (expectedStrand == ReadStrandedness::A and !isForward));
                    } else if (isForward) { // IU, ISR, OU, OSR, MU, MSR
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::A);
                    } else { // IU, ISF, OU, OSF, MU, MSF
                        return (expectedStrand == ReadStrandedness::U or
                                expectedStrand == ReadStrandedness::S);
                    }
                    break;
                default:
                    // SHOULD NOT GET HERE
                    fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
                    return false;
                    break;
            }
            // SHOULD NOT GET HERE
            fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
            return false;
        }


        // for paired-end reads
        bool compatibleHit(LibraryFormat expected, LibraryFormat observed) {
            if (observed.type != ReadType::PAIRED_END) {
                // SHOULD NOT GET HERE
                fmt::print(stderr, "WARNING: PE compatibility function called with SE read!\n");
                return false;
            }

            auto es = expected.strandedness;
            auto eo = expected.orientation;

            auto os = observed.strandedness;
            auto oo = observed.orientation;

            // If the orientations are different, they are incompatible
            if (eo != oo) {
                return false;
            } else { // In this branch, the orientations are always compatible
                return (es == ReadStrandedness::U or
                        es == os);
            }
            // SHOULD NOT GET HERE
            fmt::print(stderr, "WARNING: Could not determine strand compatibility!");
            fmt::print(stderr, "please report this.\n");
            return false;
        }


        // Determine the library type of paired-end reads
        LibraryFormat hitType(int32_t end1Start, bool end1Fwd, uint32_t len1,
                              int32_t end2Start, bool end2Fwd, uint32_t len2, bool canDovetail) {

            // If the reads come from opposite strands
            if (end1Fwd != end2Fwd) {
                // and if read 1 comes from the forward strand
                if (end1Fwd) {
                    // then if read 1 start < read 2 start ==> ISF
                    // NOTE: We can't really delineate between inward facing reads that stretch
                    // past each other and outward facing reads --- the purpose of stretch is to help
                    // make this determinateion.
                    int32_t stretch = canDovetail ? len2 : 0;
                    if (end1Start <= end2Start + stretch) {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
                    } // otherwise read 2 start < read 1 start ==> OSF
                    else {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
                    }
                }
                // and if read 2 comes from the forward strand
                if (end2Fwd) {
                    // then if read 2 start <= read 1 start ==> ISR
                    // NOTE: We can't really delineate between inward facing reads that stretch
                    // past each other and outward facing reads --- the purpose of stretch is to help
                    // make this determinateion.
                    int32_t stretch = canDovetail ? len1 : 0;
                    if (end2Start <= end1Start + stretch) {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
                    } // otherwise, read 2 start > read 1 start ==> OSR
                    else {
                        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
                    }
                }
            } else { // Otherwise, the reads come from the same strand
                if (end1Fwd) { // if it's the forward strand ==> MSF
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
                } else { // if it's the reverse strand ==> MSR
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
                }
            }
            // SHOULD NOT GET HERE
            spdlog::get("jointLog")->error("ERROR: Could not associate any known library type with read! "
                                           "Please report this bug!\n");
            std::this_thread::sleep_for(std::chrono::seconds(1));
            std::exit(-1);
            return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);
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
                ExpressionRecord(const std::string& targetIn, uint32_t lengthIn, double effLengthIn,
                        std::vector<double>& expValsIn) :
                    target(targetIn), length(lengthIn), effLength(effLengthIn), expVals(expValsIn) {}

                ExpressionRecord( ExpressionRecord&& other ) {
                    std::swap(target, other.target);
                    length = other.length;
		    effLength = other.effLength;
                    std::swap(expVals, other.expVals);
                }

                ExpressionRecord(std::vector<std::string>& inputLine) {
                    if (inputLine.size() < 3) {
                        std::string err ("Any expression line must contain at least 3 tokens");
                        throw std::invalid_argument(err);
                    } else {
                        auto it = inputLine.begin();
                        target = *it; ++it;
                        length = std::stoi(*it); ++it;
			effLength = std::stod(*it); ++it;
                        for (; it != inputLine.end(); ++it) {
                            expVals.push_back(std::stod(*it));
                        }
                    }
                }

                std::string target;
                uint32_t length;
		double effLength;
                std::vector<double> expVals;
        };

        // From : http://stackoverflow.com/questions/9435385/split-a-string-using-c11
        std::vector<std::string> split(const std::string& str, int delimiter(int) = ::isspace){
            using namespace std;
            vector<string> result;
            auto e=str.end();
            auto i=str.begin();
            while (i != e) {
                i = find_if_not(i,e, delimiter);
                if (i == e) break;
                auto j = find_if(i,e, delimiter);
                result.push_back(string(i,j));
                i = j;
            }
            return result;
        }

        /**
         * Computes (and returns) new effective lengths for the transcripts
         * based on the current abundance estimates (alphas) and the current
         * effective lengths (effLensIn).  This approach is based on the one
         * taken in Kallisto (for sequence-specific bias), and seems to work
         * well given its low computational requirements.  The handling of
         * fragment GC bias below is similar, but is not done in Kallisto.
         */
        template <typename AbundanceVecT>
          Eigen::VectorXd updateEffectiveLengths(
              SailfishOpts& sfopts,
              ReadExperiment& readExp,
              Eigen::VectorXd& effLensIn,
              AbundanceVecT& alphas) {
            using std::vector;
            double minAlpha = 1e-8;

            uint32_t gcSamp{sfopts.pdfSampFactor};
            bool gcBiasCorrect{sfopts.gcBiasCorrect};
            bool seqBiasCorrect{sfopts.biasCorrect};

            int64_t numFwdMappings = readExp.numFwd();
            int64_t numRCMappings = readExp.numRC();
            int64_t numMappings = numFwdMappings + numRCMappings;

            if (numMappings == 0) {
              sfopts.jointLog->warn("Had no fragments from which to estimate "
                  "fwd vs. rev-comp mapping rate.  Skipping "
                  "bias / gc-bias correction");
              return effLensIn;
            }

            // this should not happen for the time being
            if (gcBiasCorrect and seqBiasCorrect) {
              sfopts.jointLog->warn("enabling both sequence-specific and fragment gc bias correction "
                  "simultaneously is not yet supported.  skipping effective length "
                  "correction");
              return effLensIn;
            }

            double probFwd = static_cast<double>(numFwdMappings) / numMappings;
            double probRC = static_cast<double>(numRCMappings) / numMappings;

            // calculate read bias normalization factor -- total count in read
            // distribution.
            auto& readBias = readExp.readBias();
            int32_t K = readBias.getK();
            double readNormFactor = static_cast<double>(readBias.totalCount());

            // The *expected* biases from sequence-specific effects
            auto& transcriptKmerDist = readExp.expectedSeqBias();

            // Reset the transcript (normalized) counts
            transcriptKmerDist.clear();
            transcriptKmerDist.resize(constExprPow(4, K), 1.0);

            EmpiricalDistribution& fld = *(readExp.fragLengthDist());

            // The *expected* biases from GC effects
            auto& transcriptGCDist = readExp.expectedGCBias();
            auto& gcCounts = readExp.observedGC();
            double readGCNormFactor = 0.0;
            int32_t fldLow{0};
            int32_t fldHigh{1};

            if (gcBiasCorrect) {
              transcriptGCDist.clear();
              transcriptGCDist.resize(101, 1.0);

              bool first{false};
              bool second{false};
              for (size_t i = 0; i <= fld.maxValue(); ++i) {
                auto density = fld.cdf(i);
                if (!first and density >= 0.005) {
                  first = true;
                  fldLow = i;
                }
                if (!second and density >= 0.995) {
                  second = true;
                  fldHigh = i;
                }
              }

              for (auto& c : gcCounts) { readGCNormFactor += c; }
            }

            // Make this const so there are no shenanigans
            const auto& transcripts = readExp.transcripts();

            // The effective lengths adjusted for bias
            Eigen::VectorXd effLensOut(effLensIn.size());

            // How much to cut off
            int32_t trunc = K;

            for(size_t it=0; it < transcripts.size(); ++it) {
              auto& txp = transcripts[it];

              // First in the forward direction
              int32_t refLen = static_cast<int32_t>(txp.RefLength);
              int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

              // How much of this transcript (beginning and end) should
              // not be considered
              int32_t unprocessedLen = std::max(0, refLen - elen);

              // Skip transcripts with trivial expression or that are too
              // short.
              if (alphas[it] < minAlpha or unprocessedLen <= 0) {
                continue;
              }

              // Otherwise, proceed with the following weight.
              double contribution = (alphas[it]/effLensIn(it));

              // This transcript's sequence
              const char* tseq = txp.Sequence();

              // From the start of the transcript up until the last valid
              // kmer.
              bool firstKmer{true};
              uint32_t idx{0};

              // For each position along the transcript
              // (considering the forward direction for sequence-specific bias).
              for (int32_t i = refLen - trunc - 1; i >= 0; --i) {
                // Seq bias
                if (seqBiasCorrect) {
                  int32_t kmerStartPos = i;
                  int32_t fragStartPos = i + 2;
                  if (firstKmer) {
                    idx = indexForKmer(tseq + i, K, Direction::REVERSE_COMPLEMENT);
                    firstKmer = false;
                  } else {
                    idx = nextKmerIndex(idx, tseq[kmerStartPos], K, Direction::REVERSE_COMPLEMENT);
                  }

                  int32_t maxFragLen = refLen - fragStartPos + 1;
                  if (maxFragLen >= 0 and maxFragLen < refLen) {
                    transcriptKmerDist[idx] += probFwd * contribution * fld.cdf(maxFragLen);
                  }
                }

                // fragment GC bias
                if (gcBiasCorrect) {
                  //for (int32_t fl = fldLow; fl <= fldHigh; ++fl) {
                  double prevFLMass = fld.cdf(0);
                  for (int32_t fl = fldLow; fl <= fldHigh; fl += gcSamp) {
                    int32_t fragStart = i;
                    int32_t fragEnd = i + fl - 1; // -1 because the interval is closed on both sides
                    if (fragEnd < refLen) {
                      auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                      transcriptGCDist[gcFrac] += contribution * (fld.cdf(fl) - prevFLMass);
                      prevFLMass = fld.cdf(fl);
                    } else { break; } // no more valid positions
                  } // for each fragment  length
                } // end fragment GC bias
                } // for every position a fragment could start



                // Then in the reverse complement direction
                firstKmer = true;
                idx = 0;
                if (seqBiasCorrect) {
                  for (int32_t i = 0; i <= refLen - trunc - 1; ++i) {
                    int32_t kmerStartPos = i;
                    int32_t kmerEndPos = kmerStartPos + K - 1; // -1 because pos is *inclusive*
                    int32_t fragStartPos = kmerStartPos + 4;
                    if (firstKmer) {
                      idx = indexForKmer(tseq, K, Direction::FORWARD);
                      firstKmer = false;
                    } else {
                      idx = nextKmerIndex(idx, tseq[kmerEndPos], K, Direction::FORWARD);
                    }

                    int32_t maxFragLen = fragStartPos + 1;
                    if (maxFragLen >= 0 and maxFragLen < refLen) {
                      transcriptKmerDist[idx] += probRC * contribution * fld.cdf(maxFragLen);
                    }
                  } // end for pos in transcript
                }

              } // end for each transcript


              // Compute appropriate priors and normalization factors
              double txomeGCNormFactor = 0.0;
              double gcPrior = 0.0;
              if (gcBiasCorrect) {
                for (auto m : transcriptGCDist) { txomeGCNormFactor += m; }
                auto pmass = 101.0;
                gcPrior = ((pmass / (readGCNormFactor - pmass)) * txomeGCNormFactor) / 101.0;
              }

              double txomeNormFactor = 0.0;
              double seqPrior = 0.0;
              if (seqBiasCorrect) {
                for(auto m : transcriptKmerDist) { txomeNormFactor += m; }
                double pmass = static_cast<double>(constExprPow(4, K));
                seqPrior = ((pmass / (readNormFactor - pmass)) * txomeNormFactor) / pmass;
              }

              size_t numCorrected = 0;
              size_t numUncorrected = 0;

              // Now, compute the effective length of each transcript using
              // the k-mer biases
              for(size_t it = 0; it < transcripts.size(); ++it) {
                // Starts out as 0
                double effLength = 0.0;

                auto& txp = transcripts[it];
                // First in the forward direction, from the start of the
                // transcript up until the last valid kmer.
                int32_t refLen = static_cast<int32_t>(txp.RefLength);
                int32_t elen = static_cast<int32_t>(txp.EffectiveLength);

                // How much of this transcript (beginning and end) should
                // not be considered
                int32_t unprocessedLen = std::max(0, refLen - elen);

                Eigen::VectorXd seqFactors(refLen);
                seqFactors.setZero();
                Eigen::VectorXd gcFactors(refLen);
                gcFactors.setZero();

                if (alphas[it] >= minAlpha and unprocessedLen > 0) {
                  bool firstKmer{true};
                  uint32_t idx{0};
                  // This transcript's sequence
                  const char* tseq = txp.Sequence();

                  for (int32_t i = refLen - trunc - 1; i >= 0; --i) {
                    /** Seq-specific bias **/
                    if (seqBiasCorrect) {
                      int32_t kmerStartPos = i;
                      int32_t fragStartPos = kmerStartPos + 2;
                      if (firstKmer) {
                        idx = indexForKmer(tseq + i, K, Direction::REVERSE_COMPLEMENT);
                        firstKmer = false;
                      } else {
                        idx = nextKmerIndex(idx, tseq[kmerStartPos], K, Direction::REVERSE_COMPLEMENT);
                      }

                      int32_t maxFragLen = refLen - fragStartPos + 1;
                      if (fragStartPos >=0 and fragStartPos < refLen) {
                        seqFactors[fragStartPos] +=
                          probFwd *
                          (readBias.counts[idx]/ (transcriptKmerDist[idx] + seqPrior)) *
                          fld.cdf(maxFragLen);
                      }
                    }
                    if (gcBiasCorrect) {
                      double prevFLMass = fld.cdf(0);
                      //for (int32_t fl = fldLow; fl <= fldHigh; ++fl) {
                      for (int32_t fl = fldLow; fl <= fldHigh; fl += gcSamp) {
                        int32_t fragStart = i;
                        int32_t fragEnd = i + fl - 1; // -1 because the interval is closed on both sides
                        if (fragEnd < refLen) {
                          auto gcFrac = txp.gcFrac(fragStart, fragEnd);
                          double sampleProb = (gcCounts[gcFrac] / (gcPrior + transcriptGCDist[gcFrac])) * (fld.cdf(fl) - prevFLMass);
                          prevFLMass = fld.cdf(fl);

                          // count it in the forward orientation
                          gcFactors[fragStart] += sampleProb * probFwd;
                          // count it in the reverse compliment orientation
                          gcFactors[fragEnd] += sampleProb * probRC;
                        } else { break; } // no more valid positions
                      }
                    } // end GC bias
                    } // end checkFwd

                    // Then in the reverse complement direction
                    double seqNormFactor{0.0};
                    firstKmer = true;
                    idx = 0;
                    if (seqBiasCorrect) {
                      for (int32_t i = 0; i <= refLen - trunc - 1; ++i) {
                        int32_t kmerStartPos = i;
                        int32_t kmerEndPos = kmerStartPos + K - 1; // -1 because pos is *inclusive*
                        int32_t fragStartPos = kmerStartPos + 4;
                        if (firstKmer) {
                          idx = indexForKmer(tseq, K, Direction::FORWARD);
                          firstKmer = false;
                        } else {
                          idx = nextKmerIndex(idx, tseq[kmerEndPos], K, Direction::FORWARD);
                        }

                        int32_t maxFragLen = fragStartPos + 1;
                        if (fragStartPos >= 0 and fragStartPos < refLen) {
                          seqFactors[fragStartPos] +=
                            probRC *
                            (readBias.counts[idx]/(transcriptKmerDist[idx]+seqPrior)) *
                            fld.cdf(maxFragLen);
                        }
                      }
                      seqNormFactor = txomeNormFactor / readNormFactor;
                    }

                    if (seqBiasCorrect and gcBiasCorrect) {
                      effLength = seqFactors.dot(gcFactors);
                      effLength *= (txomeNormFactor / readNormFactor);
                      effLength *= (txomeGCNormFactor / readGCNormFactor);
                    } else if (seqBiasCorrect) {
                      effLength = seqFactors.sum();
                      effLength *= (txomeNormFactor / readNormFactor);
                    } else if (gcBiasCorrect) {
                      effLength = gcFactors.sum();
                      effLength *= (txomeGCNormFactor / readGCNormFactor);
                    }

                  } // for the processed transcript

                  if(unprocessedLen > 0.0 and effLength > unprocessedLen) {
                    ++numCorrected;
                    effLensOut(it) = effLength;
                  } else {
                    ++numUncorrected;
                    effLensOut(it) = effLensIn(it);
                  }
                }
                return effLensOut;
              }


        void aggregateEstimatesToGeneLevel(TranscriptGeneMap& tgm, boost::filesystem::path& inputPath) {

            using std::vector;
            using std::string;
            using std::ofstream;
            using std::unordered_map;
            using std::move;
            using std::cerr;
            using std::max;

            constexpr double minTPM = std::numeric_limits<double>::denorm_min();


            std::ifstream expFile(inputPath.string());

            if (!expFile.is_open()) {
                perror("Error reading file");
            }

            //====================== From GeneSum ====================
            vector<string> comments;
            unordered_map<string, vector<ExpressionRecord>> geneExps;
            string l;
            size_t ln{0};

	    bool headerLine{true};
            while (getline(expFile, l)) {
                if (++ln % 1000 == 0) {
                    cerr << "\r\rParsed " << ln << " expression lines";
                }
                auto it = find_if(l.begin(), l.end(),
                        [](char c) -> bool {return !isspace(c);});
                if (it != l.end()) {
                    if (*it == '#') {
                        comments.push_back(l);
                    } else {
		      // If this isn't the first non-comment line
		      if (!headerLine) {
			vector<string> toks = split(l);
			ExpressionRecord er(toks);
			auto gn = tgm.geneName(er.target);
			geneExps[gn].push_back(move(er));
		      } else { // treat the header line as a comment
			comments.push_back(l);
			headerLine = false;
		      }
		    }
                }
            }
            cerr << "\ndone\n";
            expFile.close();

            cerr << "Aggregating expressions to gene level . . .";
            boost::filesystem::path outputFilePath(inputPath);
            outputFilePath.replace_extension(".genes.sf");
            ofstream outFile(outputFilePath.string());

            // preserve any comments in the output
            for (auto& c : comments) {
                outFile << c << '\n';
            }

            for (auto& kv : geneExps) {
                auto& gn = kv.first;

                double geneLength = kv.second.front().length;
                double geneEffLength = kv.second.front().effLength;
                vector<double> expVals(kv.second.front().expVals.size(), 0);
                const size_t NE{expVals.size()};

                size_t tpmIdx{0};
                double totalTPM{0.0};
                for (auto& tranExp : kv.second) {
                    // expVals[0] = TPM
                    // expVals[1] = count
                    for (size_t i = 0; i < NE; ++i) { expVals[i] += tranExp.expVals[i]; }
                    totalTPM += expVals[tpmIdx];
                }

                // If this gene was expressed
                if (totalTPM > minTPM) {
                    geneLength = 0.0;
		    geneEffLength = 0.0;
                    for (auto& tranExp : kv.second) {
                        double frac = tranExp.expVals[tpmIdx] / totalTPM;
                        geneLength += tranExp.length * frac;
                        geneEffLength += tranExp.effLength * frac;
                    }
                } else {
                    geneLength = 0.0;
		    geneEffLength = 0.0;
                    double frac = 1.0 / kv.second.size();
                    for (auto& tranExp : kv.second) {
                        geneLength += tranExp.length * frac;
                        geneEffLength += tranExp.effLength * frac;
                    }
                }

                // Otherwise, if the gene wasn't expressed, the length
                // is reported as the longest transcript length.

                outFile << gn << '\t' << geneLength << '\t' << geneEffLength;
                for (size_t i = 0; i < NE; ++i) {
                    outFile << '\t' << expVals[i];
                }
                outFile << '\n';
            }

            outFile.close();
            cerr << " done\n";
            //====================== From GeneSum =====================
        }

        void generateGeneLevelEstimates(boost::filesystem::path& geneMapPath,
                boost::filesystem::path& estDir,
                std::string aggKey) {
            namespace bfs = boost::filesystem;
            std::cerr << "Computing gene-level abundance estimates\n";
            bfs::path gtfExtension(".gtf");
            auto extension = geneMapPath.extension();

            TranscriptGeneMap tranGeneMap;
            // parse the map as a GTF file
            if (extension == gtfExtension) {
                // Using libgff
                tranGeneMap = sailfish::utils::transcriptGeneMapFromGTF(geneMapPath.string(), aggKey);
            } else { // parse the map as a simple format files
                std::ifstream tgfile(geneMapPath.string());
                tranGeneMap = sailfish::utils::readTranscriptToGeneMap(tgfile);
                tgfile.close();
            }

            std::cerr << "There were " << tranGeneMap.numTranscripts() << " transcripts mapping to "
                << tranGeneMap.numGenes() << " genes\n";

            bfs::path estFilePath = estDir / "quant.sf";
            if (!bfs::exists(estFilePath)) {
                std::stringstream errstr;
                errstr << "Attempting to compute gene-level esimtates, but could not \n"
                    << "find isoform-level file " << estFilePath;
                throw std::invalid_argument(errstr.str());
            } else {
                sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, estFilePath);
            }

            /** Create a gene-level summary of the bias-corrected estimates as well if these exist **/
            /** No more of this bias correction
            if (haveBiasCorrectedFile) {
                bfs::path biasCorrectEstFilePath = estDir / "quant_bias_corrected.sf";
                if (!bfs::exists(biasCorrectEstFilePath)) {
                    std::stringstream errstr;
                    errstr << "Attempting to compute gene-level esimtates, but could not \n"
                        << "find bias-corrected isoform-level file " << biasCorrectEstFilePath;
                    throw std::invalid_argument(errstr.str());
                } else {
                    sailfish::utils::aggregateEstimatesToGeneLevel(tranGeneMap, biasCorrectEstFilePath);
                }
            }
            */
        }

        // === Explicit instantiations
        template Eigen::VectorXd updateEffectiveLengths<std::vector<tbb::atomic<double>>>(
		SailfishOpts& sfopts,
                ReadExperiment& readExp,
                Eigen::VectorXd& effLensIn,
                std::vector<tbb::atomic<double>>& alphas);

        template Eigen::VectorXd updateEffectiveLengths<std::vector<double>>(
		SailfishOpts& sfopts,
                ReadExperiment& readExp,
                Eigen::VectorXd& effLensIn,
                std::vector<double>& alphas);
    }
}

//=======
