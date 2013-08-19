#include <iostream>
#include <vector>

#include "jellyfish/parse_dna.hpp"
#include "jellyfish/mapped_file.hpp"
#include "jellyfish/parse_read.hpp"
#include "jellyfish/sequence_parser.hpp"
#include "jellyfish/dna_codes.hpp"
#include "jellyfish/compacted_hash.hpp"
#include "jellyfish/mer_counting.hpp"
#include "jellyfish/misc.hpp"

#include "HeptamerIndex.hpp"


int main(int argc, char* argv[]) {

  jellyfish::parse_read parser( argv+1, argv+2, 5000);

  jellyfish::parse_read::read_t* read;
  jellyfish::parse_read::thread stream{parser.new_thread()};
  HeptamerIndex hi;

  size_t merLen = 7;  
  uint64_t lshift{2 * (merLen - 1)};
  uint64_t masq{(1UL << (2 * merLen)) - 1};

  while ( (read = stream.next_read()) ) {
    const char         *start = read->seq_s;
    const char * const  end   = read->seq_e;

    size_t binMerId{0};
    size_t rMerId{0};
    size_t cmlen = 0;
    uint64_t kmer, rkmer;
    kmer = 0; rkmer = 0;
    // iterate over the read base-by-base
    while(start < end) {
      uint_t     c = jellyfish::dna_codes[static_cast<uint_t>(*start++)];

      switch(c) {
        case jellyfish::CODE_IGNORE: break;
        case jellyfish::CODE_COMMENT:
        std::cerr << "ERROR\n";

        // Fall through
        case jellyfish::CODE_RESET:
          cmlen = kmer = rkmer = 0;
          break;

          default:          
          kmer = ((kmer << 2) & masq) | c;
          rkmer = (rkmer >> 2) | ((0x3 - c) << lshift);
          // count if the kmer is valid in the forward and
          // reverse directions
          if(++cmlen >= merLen) {
            cmlen = merLen;
            std::cerr << "index is " << hi.index(kmer) << "\n";
          }
        }
      }

    }

}
