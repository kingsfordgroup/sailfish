Setting up Sailfish
===================

Requirements:
-------------

* A C++-11 compliant version of GCC.  Any version of [g++](gcc.gnu.org) >= 4.7 should work.

* [CMake](www.cmake.org).  Sailfish uses the CMake build system to check,
  fetch and install dependencies, and to compile and install Sailfish.  CMake
  is available for all major platforms (though Sailfish is currently
  unsupported on Windows.)

Installation:
-------------

After downloading the Sailfish source distribution and unpacking it, change into the
top-level directory:

~~~~
> cd Sailfish-0.1.0-Source
~~~~

Then, create an out-of-source build directory and change into it:

~~~~
> mkdir build
> cd build
~~~~

Sailfish makes extensive use of [Boost](www.boost.org).  We recommend
installing the most recent version (1.53) systemwide if possible. If Boost is not
installed on your system, the build process will fetch, compile and install it
locally.  However, if you already have a recent version of Boost available on
your system, it make sense to tell the build system to use that.

If you have Boost installed you can tell CMake where to look for it. Likewise, 
if you already have Intel's [Threading Building Blocks](http://threadingbuildingblocks.org/)
library installed, you can tell CMake where it is as well. The flags for CMake are as follows:

* -DFETCH_BOOST=TRUE --  If you don't have Boost installed (or have an
   older version of it), you can provide the FETCH_BOOST flag instead of the
   BOOST_ROOT variable, which will cause CMake to fetch and build Boost locally.

* -DBOOST_ROOT=<boostdir> -- Tells CMake where an existing installtion of Boost resides,
   and looks for the appropritate version in <boostdir>.  This is the top-level directory
   where Boost is installed (e.g. /opt/local).

* -DTBB_INSTALL_DIR=<tbbroot> -- Tells CMake where an existing installation of Intel's 
   TBB is installed (<tbbroot>), and looks for the apropriate headers and libraries
   there. This is the top-level directory where TBB is installed (e.g. /opt/local).

* -DCMAKE_INSTALL_PREFIX=<install_dir> -- <install_dir> is the directory to which you 
   wish Sailfish to be installed.  If you don't specify this option, it will be 
   installed locally in the top-level directory (i.e. the directory directly above "build").

Setting the appropriate flags, you can then run the CMake configure step as follows:

~~~~
> cmake [FLAGS] ..
~~~~

The above command is the cmake configuration step, which *should* complain if
anything goes wrong.  Next, you have to run the build step. Depending on what
libraries need to be fetched and installed, this could take a while
(specifically if the installation needs to install Boost).  To start the
build, just run make.

~~~~
> make
~~~~

If the build is successful, the appropriate executables and libraries should be created.
There are two points to note about the build process.  First, if the build system is 
downloading and compiling boost, you may see a large number of warnings during compilation;
these are normal.  Second, note that CMake has colored output by default, and the steps which
create or link libraries are printed in red.  This is the color chosen by CMake for linking 
messages, and does not denote an error in the build process. 

Finally, after everything is built, the libraries and executable can be installed with:

~~~~
> make install
~~~~

To ensure that Sailfish has access to the appropriate libraries you should ensure
that the PATH variabile contains \<install_dir\>/bin, and that LD_LIBRARY_PATH 
(or DYLD_FALLBACK_LIBRARY_PATH on OSX) contains \<install_dir\>/lib.

After the paths are set, you can test the installation by running

~~~~
> make test
~~~~

This should run a simple test and tell you if it succeeded or not.

Running Sailfish
================

Sailfish is a reference based isoform quantification tool.  Thus, it requires a reference
transcriptome to run.  All you need to run Sailfish is a fasta file containing your reference 
transcripts and a (set of) fasta/fastq file(s) containing your reads.  Sailfish runs in two
phases; indexing and quantification.  The indexing step is independent of the reads, and only
need to be run one for a particular set of reference transcripts and choice of k (the kmer size).
The quantification step, obviously, is specific to the set of RNA-seq reads and is thus run more
frequently.

Indexing
--------

To generate the Sailfish index for your reference set of transcripts, you should run the following
command:

~~~~
> sailfish index -t <ref_transcripts> -o <out_dir> -k <kmer_len>
~~~~

This will build a Sailfish index for kmers of length \<kmer_len\> for the
reference transcripts  provided in the file \<ref_transcripts\> and place the
index under the directory \<out_dir\>.  There  are additional options that can
be passed to the Sailfish indexer (e.g. the number of threads to use).  These
can be seen by executing the command "Sailfish index -h".

Quantification
--------------

Now that you have generated the Sailfish index (say that it's the directory \<index_dir\> --- this
corresponds to the \<out_dir\> argument provided in the previous step), you can quantify the transcript
expression for a given set of reads.  To perform the quantification, you run a command like the following:

~~~~
> sailfish quant -i <index_dir> --reads <reads_1> <reads_2> . . . <reads_n> -o <quant_dir>
~~~~

Where \<index_dir\> is, as described above, the location of the sailfish index, \<reads_{1 . . . n}\>
are the files containing the reads, and \<quant_dir\> is the directory where the output should be written.
Just like the indexing step, additional options are available, and can be viewed by running
"Sailfish quant -h".

When the quantification step is finished, the directory \<quant_dir\> will conatin a file named
"quant.sf".  This file contains the result of the Sailfish quantification step.  This file contains a
number of columns (which are listed in the last of the header lines beginning with '#').  Specifically,
the columns are (1) Transcript ID, (2) Transcript Length, (3) Transcripts per Million (TPM), 
and (4) Reads Per Kilobase per Million mapped reads (RPKM).  The first two columns are self-explanatory,
while the last two are measures of transcript abundance.  The Transcripts per Million quantification
number is computed as described in [1], and is meant as an estimate of the number of transcripts, per
million observed transcripts, originating from each isoform.  Its benefit over the RPKM measure is
that it is independent of the mean expressed transcript length (i.e. if the mean expressed transcript
length varies between samples, for example, this alone can affect differential analysis based on 
the RPKM.)  The RPKM is a classic measure of relative transcript abundance, and is an estimate of the
number of reads per kilobase of transcript (per million mapped reads) originating from each transcript.

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

References
==========

[1] Li, Bo, et al. "RNA-Seq gene expression estimation with read mapping uncertainty." 
    Bioinformatics 26.4 (2010): 493-500.

















