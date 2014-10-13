#include <boost/thread/thread.hpp>

#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <atomic>
#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <cassert>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>


#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include "jellyfish/config.h"
#include "jellyfish/err.hpp"
#include "jellyfish/misc.hpp"
#include "jellyfish/jellyfish.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"


#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

#include "tbb/parallel_for_each.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_sort.h"

#if HAVE_LOGGER
#include "g2logworker.h"
#include "g2log.h"
#endif

#include "cmph.h"
#include "CountDBNew.hpp"
#include "LookUpTableUtils.hpp"
#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "PerfectHashIndex.hpp"

using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

typedef jellyfish::stream_manager<char**>                stream_manager;
typedef jellyfish::whole_sequence_parser<stream_manager> sequence_parser;

extern "C" {
int bwa_index(int argc, char* argv[]);
}

int salmonIndex(int argc, char* argv[]) {

    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<string>()->required(), "Transcript fasta file.")
    ("index,i", po::value<string>()->required(), "Salmon index.")
    ;

    po::variables_map vm;
    int ret = 1;
    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Index
==========
Augments an existing Sailfish index [index] with a
Salmon index if it exists, or creates a new index
[index] if one does not.
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        string transcriptFiles = vm["transcripts"].as<string>();
        bfs::path indexDirectory(vm["index"].as<string>());

        if (!bfs::exists(indexDirectory)) {
            std::cerr << "index [" << indexDirectory << "] did not previously exist "
                      << " . . . creating it\n";
            bfs::create_directory(indexDirectory);
        }

        bfs::path outputPrefix = indexDirectory / "bwaidx";

        std::vector<char*> bwaArgVec{ "index", "-p",
                                    const_cast<char*>(outputPrefix.string().c_str()),
                                    const_cast<char*>(transcriptFiles.c_str()) };

        char* bwaArgv[] = { bwaArgVec[0], bwaArgVec[1],
                            bwaArgVec[2], bwaArgVec[3] };
        int bwaArgc = 4;

        ret = bwa_index(bwaArgc, bwaArgv);

    std::cerr << "done\n";

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " index was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }
    return ret;
}

/*
int salmonIndex(int argc, char* argv[]) {

    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("transcripts,t", po::value<string>(), "Transcript fasta file.")
    ("index,i", po::value<string>(), "sailfish index.")
    ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use concurrently.")
    ;

    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(generic).run(), vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Index
==========
Augments an existing Sailfish index [index] with a
Salmon index.
)";
            std::cout << hstring << std::endl;
            std::cout << generic << std::endl;
            std::exit(1);
        }
        po::notify(vm);

        string transcriptFiles = vm["transcripts"].as<string>();
        uint32_t nbThreads = vm["threads"].as<uint32_t>();
        bfs::path indexDirectory(vm["index"].as<string>());

        const int concurrentFile = 1;
        const int maxReadGroup = 10;

        char* transcriptFileArray[] = { const_cast<char*>(transcriptFiles.c_str()) };
        stream_manager streams(transcriptFileArray, transcriptFileArray + 1, concurrentFile);

        sequence_parser parser(4 * nbThreads, maxReadGroup, concurrentFile, streams);

        bfs::path sfIndexPath = indexDirectory / "transcriptome.sfi";
        std::string sfTrascriptIndexFile(sfIndexPath.string());
        std::cerr << "reading index . . . ";
        auto phi = PerfectHashIndex::fromFile(sfTrascriptIndexFile);
        std::cerr << "done\n";
        std::cerr << "index contained " << phi.numKeys() << " kmers\n";

        std::unordered_map<std::string, uint64_t> transcriptIndexMap;

         bfs::path tlutPath = indexDirectory / "transcriptome.tlut";
         // Get transcript lengths
         std::ifstream ifile(tlutPath.string(), std::ios::binary);
         size_t numRecords {0};
         ifile.read(reinterpret_cast<char *>(&numRecords), sizeof(numRecords));

         std::cerr << "Transcript LUT contained " << numRecords << " records\n";
         for (auto i : boost::irange(size_t(0), numRecords)) {
             auto ti = LUTTools::readTranscriptInfo(ifile);
             // copy over the length, then we're done.
             transcriptIndexMap[ti->name] = ti->transcriptID;
         }
         ifile.close();
         // --- done ---

         size_t nkeys = phi.numKeys();
         size_t merLen = phi.kmerLength();
         std::cerr << "kmer length = " << merLen << "\n";

        // Read in the Jellyfish hash of the transcripts
         bfs::path thashFile = indexDirectory / "jf.counts";
         std::ifstream transcriptDB(thashFile.string());
         if (!transcriptDB.good()) {
#if HAVE_LOGGER
             LOG(FATAL) << "Couldn't open the Jellyfish hash [" << thashFile << "] quitting\n";
#endif
             std::exit(-1);
         }
         jellyfish::file_header header;
         header.read(transcriptDB);

         std::cerr << "transcript hash size is " << header.size() << "\n";
         //size_t nkeys = header.size();
         // Since JF2, key_len() is in terms of bits so merLen = keyLen / 2;
         //size_t merLen = header.key_len() / 2;
         jellyfish::mer_dna::k(merLen);

         std::cerr << "header.key_len = " << merLen << "\n";
         std::vector<std::tuple<uint64_t, uint32_t>> khits;
         size_t totHits{0};

         // Read in the jellyfish hash to determine how many hits there are
         // for each k-mer
         if (!header.format().compare(binary_dumper::format)) {
             binary_reader reader(transcriptDB, &header);
             while ( reader.next() ) {
                 auto key = reader.key().get_bits(0, 2*merLen);
                 auto id = phi.index(key);
                 auto numHits = reader.val();
                 khits.emplace_back(std::forward_as_tuple(id, numHits));
                 totHits += numHits;
             }

         } else if (!header.format().compare(text_dumper::format)) {
             text_reader reader(transcriptDB, &header);
             while ( reader.next() ) {
                 auto key = reader.key().get_bits(0, 2*merLen);
                 auto id = phi.index(key);
                 auto numHits = reader.val();
                 khits.emplace_back(std::forward_as_tuple(id, numHits));
                 totHits += numHits;
             }
         } else {
#if HAVE_LOGGER
             LOG(FATAL) << "Unknown Jellyfish hash format. quitting\n";
#endif
             std::exit(-1);
         }

         // Sort the k-mer's by their perfect hash ID
         std::cerr << "sorting \n";
         tbb::parallel_sort(khits.begin(), khits.end(),
                 [](const std::tuple<uint64_t, uint32_t>& t1, const std::tuple<uint64_t, uint32_t>& t2) -> bool {
                 return std::get<0>(t1) < std::get<0>(t2);
                 });
         std::cerr << "done";

         std::vector<uint64_t> offsets(nkeys + 1, 0);
         std::vector<std::atomic<uint32_t>> curNum(nkeys);
         size_t curOffset{0};
         for (auto& idOff : khits) {
             auto id = std::get<0>(idOff);
             auto count = std::get<1>(idOff);
             offsets[id] = curOffset;
             curOffset += count;
         }
         offsets[nkeys] = curOffset;

         std::cerr << "curOffset = " << curOffset << ", totalHits = " << totHits << "\n";
         std::vector<uint64_t> kmerMap(totHits, 0);

         auto INVALID = phi.INVALID;

         size_t effThreads = nbThreads > 1 ? (nbThreads - 1) : 1;

         std::atomic<uint32_t> tid{0};

         std::vector<std::thread> parseThreads;
         for (size_t threadNum = 0; threadNum < effThreads; ++threadNum) {

             parseThreads.push_back( std::thread([&]() -> void {
               while(true) {
                sequence_parser::job j(parser); // Get a job from the parser: a bunch of read (at most max_read_group)
                if(j.is_empty()) break;          // If got nothing, quit

                 my_mer kmer;

                 for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got
                     const char* start     = j->data[i].seq.c_str();
                     uint32_t readLen      = j->data[i].seq.size();
                     const char* const end = start + readLen;

                     // The transcript name
                     std::string fullHeader(j->data[i].header);
                     std::string tname = fullHeader.substr(0, fullHeader.find(' '));
                     uint32_t cmlen{0};

                     uint32_t rbase{0};
                     uint32_t tnum = transcriptIndexMap[tname];
                     // iterate over the read base-by-base
                     while(start < end) {
                         ++rbase;
                         char base = *start; ++start;
                         auto c = jellyfish::mer_dna::code(base);
                         kmer.shift_left(c);

                         switch(c) {
                             case jellyfish::mer_dna::CODE_IGNORE:
                                 break;
                             case jellyfish::mer_dna::CODE_COMMENT:
                                 std::cerr << "ERROR: unexpected character " << base << " in read!\n";
                                 // Fall through
                             case jellyfish::mer_dna::CODE_RESET:
                                 cmlen = 0;
                                 kmer.polyA();
                                 break;

                             default:
                                 if(++cmlen >= merLen) {
                                     cmlen = merLen;

                                     // Get the index of the current k-mer
                                     auto merID = phi.index(kmer.get_bits(0, 2*merLen));

                                     // If it's a valid k-mer (shouldn't it always be so here?)
                                     if (merID != INVALID) {

                                         // We're in transcript tnum, and the k-mer
                                         // begins at position rbase - merLen
                                         auto enc = sailfish::utils::encode(tnum, rbase - merLen);

                                         // The offset of the hit list of this k-mer
                                         auto koffset = offsets[merID];
                                         // The index of this hit in the hit list
                                         auto curOffset = curNum[merID]++;
                                         // Put this hit in the list
                                         kmerMap[koffset + curOffset] = enc;
                                     } // merID != INALID
                                 } // ++cmlen >= merLen
                         } // end switch(c)
                     } // while (start < end)

                     // Report our progress every so often
                     ++tid;
                     if (tid % 1000 == 0) {
                         std::cerr << "\r\rprocessed " << tid << " transcripts";
                     }
                 }
            } // while(true)
        }));
    } // creating parsing threads

    // wait for parsing threads to finish
    for (auto& t : parseThreads) { t.join(); }
    std::cerr << "done\n";

    tbb::parallel_for(tbb::blocked_range<size_t>(0, nkeys),
                    [&] (tbb::blocked_range<size_t>& krange) -> void {

                    for (auto kidx = krange.begin(); kidx != krange.end(); ++kidx) {
                        auto count = offsets[kidx+1] - offsets[kidx];
                        std::sort(kmerMap.begin() + offsets[kidx], kmerMap.begin() + offsets[kidx] + count);
                    }
    });

    bfs::path lookupPath = indexDirectory / "fullLookup.kmap";

    std::cerr << "writing salmon index to " << lookupPath << " . . . ";
    std::ofstream binstream(lookupPath.string(), std::ios::binary);
    cereal::BinaryOutputArchive archive(binstream);
    archive(offsets, kmerMap);
    binstream.close();
    std::cerr << "done\n";

    } catch (po::error &e) {
        std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " index was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " index --help\nExiting.\n";
    }
}
*/
