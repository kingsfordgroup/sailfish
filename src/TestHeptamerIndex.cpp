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
