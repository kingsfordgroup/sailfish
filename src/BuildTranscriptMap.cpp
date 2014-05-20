
#include <boost/thread/thread.hpp>
#include <boost/lockfree/queue.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <string>
#include <memory>
#include <functional>
#include <unordered_map>
#include <mutex>
#include <thread>
#include <chrono>
#include <iomanip>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/range/irange.hpp>
#include <boost/filesystem.hpp>

// Jellyfish includes
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "SailfishUtils.hpp"
#include "GenomicFeature.hpp"
#include "CountDBNew.hpp"
#include "SailfishConfig.hpp"
#include "VersionChecker.hpp"

bool parseTranscripts(boost::filesystem::path& tpath, PerfectHashIndex& phi, boost::filesystem::path& opath) {
    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;
    using stream_manager = jellyfish::stream_manager<char**>;
    using sequence_parser = jellyfish::whole_sequence_parser<stream_manager>;

    using BinMer = uint64_t;
    auto INVALID = phi.INVALID;
    auto merLen = phi.kmerLength();

    // set the appropriate k-mer length
    my_mer::k(merLen);

    std::ofstream ofile(opath.string(), std::ios::out | std::ios::binary);

    int l;
    uint64_t seqnum{0};

    // Create a jellyfish parser
    char** fnames = new char*[1];// fnames[1];
    fnames[0] = const_cast<char*>(tpath.c_str());

    const int concurrentFile{1};
    stream_manager streams(fnames, fnames + 1, concurrentFile);

    const uint32_t nworking{1};
    size_t maxReadGroupSize{100};
    sequence_parser parser(4*nworking, maxReadGroupSize, concurrentFile, streams);



    int fp = open(tpath.string().c_str(), O_RDONLY);
    BinMer masq{(1UL << (2 * merLen)) - 1};

    //seq = kseq_init(fp);

    std::vector<BinMer> klist;
    while(true) {
        sequence_parser::job j(parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;          // If got nothing, quit

        my_mer kmer;

        for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read we got

            char* nameStart = const_cast<char*>(j->data[i].header.c_str());
            uint32_t nameLen = j->data[i].header.size();
            char* nameEnd = nameStart + nameLen;

            char* start     = const_cast<char*>(j->data[i].seq.c_str());
            uint32_t readLen      = j->data[i].seq.size();
            char* const end = start + readLen;

            /*
    while ((l = kseq_read(seq)) >= 0) {

        char* nameStart = seq->name.s;
        uint32_t nameLen = seq->name.l;
        char* nameEnd = seq->name.s + seq->name.l;
        char* start     = seq->seq.s;
        uint32_t readLen      = seq->seq.l;
        char* end = seq->seq.s + readLen;
            */

        // reset
        uint64_t cmlen{0};
        uint64_t kmer{0};
        uint64_t numKmers{0};

        if (readLen > klist.size()) {
            klist.reserve(readLen);
        }

        size_t binMerId{0};
        // iterate over the read base-by-base
        while(start < end) {
            char pc = *start++;
            if (pc == ' ' or pc == '\n') { continue; }
            auto c = jellyfish::mer_dna::code(pc);
            //kmer.shift_left(c);

            switch(c) {
                case jellyfish::mer_dna::CODE_IGNORE: break;
                case jellyfish::mer_dna::CODE_COMMENT:
                std::cerr << "ERROR: unexpected character " << c << " in read!\n";
                // Fall through
                case jellyfish::mer_dna::CODE_RESET:
                cmlen = kmer = 0;
                break;
            default:
                kmer = ((kmer << 2) & masq) | c;
                if(++cmlen >= merLen) {
                    cmlen = merLen;

                    binMerId = phi.index(kmer);
                    if (binMerId == INVALID) {
                        std::cerr << "last base = " << static_cast<char>(c) << "\n";
                        throw std::invalid_argument("found invalid kmer in transcript");
                    }
                    klist[numKmers] = binMerId;
                    ++numKmers;
                } // if k-mer is long enough
            } // end switch
        } // all k-mers in this sequence


        ofile.write(reinterpret_cast<char*>(&nameLen), sizeof(nameLen));
        ofile.write(reinterpret_cast<char*>(nameStart), nameLen * sizeof(char));
        ofile.write(reinterpret_cast<char*>(&numKmers), sizeof(numKmers));
        ofile.write(reinterpret_cast<char*>(&klist[0]), numKmers * sizeof(klist[0]));

        ++seqnum;
        if (seqnum % 10000 == 0) {
            std::cerr << "processed: " << seqnum << " transcripts\r\r";
        }
    }
    } // end while(true)
    std::cerr << "\n";
    ofile.close();
    //kseq_destroy(seq);
    close(fp);
}

int main(int argc, char* argv[] ) {
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  string cmdString = "tmap";

  try{

   uint32_t maxThreads = std::thread::hardware_concurrency();

    po::options_description generic("Command Line Options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
    ;

    po::options_description config("Configuration");
    config.add_options()
      ("index,i", po::value<string>(), "sailfish index prefix (without .sfi/.sfc)")
      ("transcripts,t", po::value<string>(), "file containing the transcripts")
      ("output,o", po::value<string>(), "output file")
        ("threads,p", po::value<uint32_t>()->default_value(maxThreads), "The number of threads to use when counting kmers")
      ;

    po::options_description programOptions("combined");
    programOptions.add(generic).add(config);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(programOptions).run(), vm);

    //bool poisson = ( vm.count("poisson") ) ? true : false;

    if ( vm.count("version") ) {
      std::cout << "version : " << sailfish::version <<"\n";
      std::exit(0);
    }

    if ( vm.count("help") ){
      std::cout << "Sailfish\n";
      std::cout << programOptions << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    uint32_t numThreads = vm["threads"].as<uint32_t>();
    //tbb::task_scheduler_init init(numThreads);

    string transcriptFileName = vm["transcripts"].as<string>();
    bfs::path transcriptFile(transcriptFileName);

    string sfIndexBase = vm["index"].as<string>();
    bfs::path sfIndexBasePath(vm["index"].as<string>());

    string sfIndexFile = sfIndexBase+".sfi";
    string sfTrascriptCountFile = sfIndexBase+".sfc";
    bfs::path outputFilePath = bfs::path(vm["output"].as<string>());

    auto kmerEquivClassFname = bfs::path(sfIndexBase);
    kmerEquivClassFname = kmerEquivClassFname.parent_path();
    kmerEquivClassFname /= "kmerEquivClasses.bin";

    std::cerr << "Reading transcript index from [" << sfIndexFile << "] . . .";
    auto sfIndex = PerfectHashIndex::fromFile( sfIndexFile );
    std::cerr << "done\n";
    parseTranscripts(transcriptFile, sfIndex, outputFilePath);

  } catch (po::error &e){
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::cerr << argv[0] << " " << cmdString << " was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0] << " " << cmdString << " --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}

