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

// holding 2-mers as a uint64_t is a waste of space,
// but using Jellyfish makes life so much easier, so
// we'll live with it for now.
using Kmer = uint64_t;
using Sailfish::TranscriptFeatures;
namespace bfs = boost::filesystem;


int computeBiasFeatures(
    std::vector<std::string>& transcriptFiles,
    bfs::path outFilePath,
    size_t numThreads) {

    using std::string;

	std::vector<std::string> alphabet{{"AA", "AC", "AG", "AT",
                        "CA", "CC", "CG", "CT",
                        "GA", "GC", "GG", "GT",
	                    "TA", "TC", "TG", "TT"}};

	std::vector<Kmer> diNucleotides(16);
	for (auto i : boost::irange(size_t{0}, diNucleotides.size())) {
		diNucleotides[i] = jellyfish::parse_dna::mer_string_to_binary(alphabet[i].c_str(), 2);
	}


        std::vector<std::string> readFiles = transcriptFiles;
        for( auto rf : readFiles ) {
            std::cerr << "readFile: " << rf << ", ";
        }
        std::cerr << "\n";

        char** fnames = new char*[readFiles.size()];
        size_t z{0};
        size_t numFnames{0};
        for ( auto& s : readFiles ){
            // Ugly, yes?  But this is not as ugly as the alternatives.
            // The char*'s contained in fname are owned by the readFiles
            // vector and need not be manually freed.
            fnames[numFnames] = const_cast<char*>(s.c_str());
            ++numFnames;
        }

        // Create a jellyfish parser
        jellyfish::parse_read parser( fnames, fnames+numFnames, 5000);


	    size_t merLen = 2;
        Kmer lshift(2 * (merLen - 1));
        Kmer masq((1UL << (2 * merLen)) - 1);
        std::atomic<size_t> readNum{0};


        size_t numActors = numThreads;
        size_t numComplete = 0;
	    std::vector<std::thread> threads;
        auto tstart = std::chrono::steady_clock::now();

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
        	}
        );

        for (auto i : boost::irange(size_t{0}, numActors)) {
	    threads.push_back(std::thread(
	        [&featQueue, &numComplete, &parser, &readNum, &tstart, lshift, masq, merLen, numActors]() -> void {

                jellyfish::parse_read::read_t* read;
                jellyfish::parse_read::thread stream = parser.new_thread();
                size_t cmlen, kmer, numKmers;
                while ( (read = stream.next_read()) ) {
                    ++readNum; //++locallyProcessedReads;
                    if (readNum % 1000 == 0) {
                        auto tend = std::chrono::steady_clock::now();
                        auto sec = std::chrono::duration_cast<std::chrono::seconds>(tend-tstart);
                        auto nsec = sec.count();
                        auto rate = (nsec > 0) ? readNum / sec.count() : 0;
                        std::cerr << "processed " << readNum << " transcripts (" << rate << ") transcripts/s\r\r";
                    }


                    // we iterate over the entire read
                    const char         *start = read->seq_s;
                    const char * const  end   = read->seq_e;

					TranscriptFeatures tfeat{};

                    // reset all of the counts
                    numKmers = 0;
                    cmlen = kmer = 0;

                    uint32_t readLen = std::distance(start, end);
                    // The transcript name
                    std::string fullHeader(read->header, read->hlen);
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

                } // end reads
            } // end lambda
            ));

		} // actor loop

		for (auto& t : threads) { t.join(); ++numComplete; }
		std::cerr << "\n";
		outputThread.join();


}
