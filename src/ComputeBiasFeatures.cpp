
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

#include <iostream>
#include <vector>
#include <array>
#include <atomic>
#include <thread>

#include "jellyfish/parse_dna.hpp"
#include "jellyfish/mapped_file.hpp"
#include "jellyfish/parse_read.hpp"
#include "jellyfish/sequence_parser.hpp"
#include "jellyfish/dna_codes.hpp"
#include "jellyfish/compacted_hash.hpp"
#include "jellyfish/mer_counting.hpp"
#include "jellyfish/misc.hpp"

#include "tbb/concurrent_queue.h"

#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "CommonTypes.hpp"
#include "ReadProducer.hpp"
#include "StreamingSequenceParser.hpp"

// holding 2-mers as a uint64_t is a waste of space,
// but using Jellyfish makes life so much easier, so
// we'll live with it for now.
using Kmer = uint64_t;
using Sailfish::TranscriptFeatures;
namespace bfs = boost::filesystem;

template <typename ParserT>
bool computeBiasFeaturesHelper(ParserT& parser,
                               tbb::concurrent_bounded_queue<TranscriptFeatures>& featQueue,
                               size_t& numComplete, size_t numThreads) {
    size_t merLen = 2;
    Kmer lshift(2 * (merLen - 1));
    Kmer masq((1UL << (2 * merLen)) - 1);
    std::atomic<size_t> readNum{0};

    size_t numActors = numThreads;
    std::vector<std::thread> threads;
    auto tstart = std::chrono::steady_clock::now();

    for (auto i : boost::irange(size_t{0}, numActors)) {
        threads.push_back(std::thread(
	        [&featQueue, &numComplete, &parser, &readNum, &tstart, lshift, masq, merLen, numActors]() -> void {

                ReadProducer<ParserT> producer(parser);

                ReadSeq* s;
                size_t cmlen, kmer, numKmers;
                while (producer.nextRead(s)) {
                    ++readNum; 
                    if (readNum % 1000 == 0) {
                        auto tend = std::chrono::steady_clock::now();
                        auto sec = std::chrono::duration_cast<std::chrono::seconds>(tend-tstart);
                        auto nsec = sec.count();
                        auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                        std::cerr << "processed " << readNum << " transcripts (" << rate << ") transcripts/s\r\r";
                    }

                    // we iterate over the entire read
                    const char* start     = s->seq;
                    uint32_t readLen      = s->len;
                    const char* const end = s->seq + readLen;

                    TranscriptFeatures tfeat{};

                    // reset all of the counts
                    numKmers = 0;
                    cmlen = kmer = 0;

                    // the maximum number of kmers we'd have to store
                    uint32_t maxNumKmers = (readLen >= merLen) ? readLen - merLen + 1 : 0;
                    if (maxNumKmers == 0) { featQueue.push(tfeat); continue; }
                    
                    // The transcript name
                    std::string fullHeader(s->name, s->nlen);
                    tfeat.name = fullHeader.substr(0, fullHeader.find(' '));
                    tfeat.length = readLen;
                    auto nfact = 1.0 / readLen;

                    // iterate over the read base-by-base
                    while(start < end) {
                    	    char base = *start;
                            uint_t     c = jellyfish::dna_codes[static_cast<uint_t>(*start++)];

                            // ***** Potentially consider quality values in the future **** /
                            // const char q = *start++;
                            // if(q < q_thresh)
                            //   c = CODE_RESET;

                            switch(c) {
                                case jellyfish::CODE_IGNORE: break;
                                case jellyfish::CODE_COMMENT:
                                  std::cerr << "ERROR\n";
                                  //report_bad_input(*(start-1));
                                // Fall through
                                case jellyfish::CODE_RESET:
                                  cmlen = kmer = 0;
                                  break;

                                default:
                                  // form the new kmer
                                  kmer = ((kmer << 2) & masq) | c;
                                  if (++cmlen >= merLen) {
                                      tfeat.diNucleotides[kmer]++;
                                      if (base == 'G' or base == 'C') { tfeat.gcContent += nfact; }
                                  }

                            } // end switch

                        } // end while

                    char lastBase = *(end - 1);
                    if (lastBase == 'G' or lastBase == 'C') { tfeat.gcContent += nfact; }
                    featQueue.push(tfeat);

                    producer.finishedWithRead(s);

                } // end reads
            } // end lambda
            ));

        } // actor loop

        for (auto& t : threads) { t.join(); ++numComplete; }
        return true;
}

int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    bfs::path outFilePath,
    bool useStreamingParser,
    size_t numThreads) {

    using std::string;
    using std::vector;
    using std::cerr;

    size_t numActors = numThreads;
    size_t numComplete = 0;
    tbb::concurrent_bounded_queue<TranscriptFeatures> featQueue;

    std::ofstream ofile(outFilePath.string());

    auto outputThread = std::thread(
         [&ofile, &numComplete, &featQueue, numActors]() -> void {
             TranscriptFeatures tf{};
             while( numComplete < numActors or !featQueue.empty() ) {
                 while(featQueue.try_pop(tf)) {
                     ofile << tf.name << '\t';
                     ofile << tf.length << '\t';
                     ofile << tf.gcContent << '\t';
                     for (auto i : boost::irange(size_t{0}, tf.diNucleotides.size())) {
                         ofile << tf.diNucleotides[i];
                         char end = (i == tf.diNucleotides.size() - 1) ? '\n' : '\t';
                         ofile << end;
                     }
                 }
                 boost::this_thread::sleep_for(boost::chrono::milliseconds(100));

             }
             ofile.close();
         });

    bool tryJellyfish = !useStreamingParser;

    std::vector<std::string> readFiles = transcriptFiles;
    for( auto rf : readFiles ) {
        std::cerr << "readFile: " << rf << ", ";
    }
    std::cerr << "\n";

    for (auto& readFile : readFiles) {
        std::cerr << "file " << readFile << ": \n";

        namespace bfs = boost::filesystem;
        bfs::path filePath(readFile);

        // If this is a regular file, then use the Jellyfish parser
        if (tryJellyfish and bfs::is_regular_file(filePath)) {

            char** fnames = new char*[1];// fnames[1];
            fnames[0] = const_cast<char*>(readFile.c_str());

            jellyfish::parse_read parser(fnames, fnames+1, 5000);

            computeBiasFeaturesHelper<jellyfish::parse_read>(
                                                             parser, featQueue, numComplete, numActors);

        } else { // If this is a named pipe, then use the kseq-based parser
            vector<bfs::path> paths{readFile};
            StreamingReadParser parser(paths);
            parser.start();
            computeBiasFeaturesHelper<StreamingReadParser>(
                                                           parser, featQueue, numComplete, numActors);
        }
    }

    std::cerr << "\n";
    outputThread.join();
}
