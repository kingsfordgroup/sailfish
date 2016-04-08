**Note**: If you're using a version of Sailfish prior to v0.7.0, please consider upgrading to the latest version as soon as is convenient.  It contains **substantial** improvements over v0.6.3 in terms of both accuracy and speed.

**Question**: I've read the [kallisto](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3519.html) paper, and it looks like kallisto is much more accurate (and faster) than the Sailfish software; is this true?

**Answer**: No. The benchmarks in the kallisto paper are made against a [**very old**](https://github.com/pachterlab/kallisto_paper_analysis/blob/nbt/config.py#L44) version of Sailfish (from March of 2014), which has now been deprecated for quite some time.  Based on the success we encountered by replacing alignment with mapping (instead of k-mer counting) in [Salmon](https://COMBINE-lab.github.io/salmon), we decided to update Sailfish to make use of this improved mapping information.  As mentioned in the note above, this results in substantial improvements in both Sailfish's accuracy and its speed.  For a comparison of quantification accuracy that considers reasonably recent versions of kallisto and Sailfish, please refer to the RapMap [pre-print](http://biorxiv.org/content/early/2015/10/28/029652) (the most-recent pre-print is available [here](http://biorxiv.org/content/early/2016/01/16/029652)).  You can also see a re-analysis of Sailfish's performance on some of the kallisto simulations (using a more recent version of Sailfish) in the [QuantAnalysis repo](https://github.com/COMBINE-lab/QuantAnalysis).

If you're interested in Sailfish, you might also want to take a look
at our new software, [Salmon](https://COMBINE-lab.github.io/salmon).

The documentation for Sailfish is being migrated to [ReadTheDocs](http://readthedocs.org).
To see [the latest documentation there](http://sailfish.readthedocs.org).

Documentation, Build Status and Gitter Channel:
===============================================

[![Documentation Status](https://readthedocs.org/projects/sailfish/badge/?version=master)](http://sailfish.readthedocs.org)
[![Build Status](https://travis-ci.org/kingsfordgroup/sailfish.svg?branch=master)](https://travis-ci.org/kingsfordgroup/sailfish)
[![Join the chat at https://gitter.im/kingsfordgroup/ailfish](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/kingsfordgroup/sailfish?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Requirements:
-------------

To build Sailfish from source, the following are required:

* A C++-11 compliant version of GCC.  Any version of [g++](gcc.gnu.org) >= 4.7
  should work.

* [CMake](www.cmake.org).  Sailfish uses the CMake build system to check,
  fetch and install dependencies, and to compile and install Sailfish.  CMake
  is available for all major platforms (though Sailfish is currently
  unsupported on Windows.)

Installation and Usage:
-----------------------

For information about building, installing or using Sailfish, please refer
to [the documentation](http://sailfish.readthedocs.org/).

License
=======

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see \<http://www.gnu.org/licenses/\>.
