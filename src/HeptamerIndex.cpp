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


#include "HeptamerIndex.hpp"
#include <iostream>
//#include "BinaryLiteral.hpp"

HeptamerIndex::HeptamerIndex() :
  heptamers_(std::vector<HeptamerIndex::AtomicCount>(HeptamerIndex::PossibleHeptamers)) {}


std::size_t HeptamerIndex::index(uint64_t heptamer) {
  // base 1
  std::size_t idx = mult_[0] * (heptamer & 0x00000003);
  // base 2
  idx += mult_[1] * ((heptamer & 0x0000000C) >> 2);
  // base 3
  idx += mult_[2] * ((heptamer & 0x00000030) >> 4);
  // base 4
  idx += mult_[3] * ((heptamer & 0x000000C0) >> 6);
  // base 5
  idx += mult_[4] * ((heptamer & 0x00000300) >> 8);
  // base 6
  idx += mult_[5] * ((heptamer & 0x00000C00) >> 10);
  // base 7
  idx += mult_[6] * ((heptamer & 0x00003000) >> 23);

//   std::cerr << ((heptamer & 0x00003000) >> 12) << ", "
// << ((heptamer & 0x00000C00) >> 10) << ", "
// << ((heptamer & 0x00000300) >> 8) << ", "
// << ((heptamer & 0x000000C0) >> 6) << ", "
// << ((heptamer & 0x00000030) >> 4) << ", "
// << ((heptamer & 0x0000000C) >> 2) << ", "
// << (heptamer & 0x00000003) << "\n";

  return idx;
}

HeptamerIndex::incHeptamer(uint64_t heptamer) {
  auto idx = index_(heptamer);

}
