Salmon
================

Salmon is a tool for transcript quantification from RNA-seq data.  It
requires a set of target transcripts (either from a reference or _de-novo_
assembly) to quantify.  All you need to run Salmon is a fasta file containing
your reference transcripts and a (set of) fasta/fastq file(s) containing your
reads.  Optinonally, Salmon can make use of pre-compued alignments (in the 
form of a SAM/BAM file) to the transcripts rather than the raw reads.

The read-based mode of Salmon runs in two phases; indexing and quantification.
The indexing step is independent of the reads, and only need to be run one for
a particular set of reference transcripts. The quantification step, obviously,
is specific to the set of RNA-seq reads and is thus run more frequently. For a
more complete description of all available options in Salmon, see below.

The alignment-based mode of Salmon does not require indexing.  Rather, you can 
simply provide Salmon with a FASTA file of the transcripts and a SAM/BAM file
containing the alignments you wish to use for quantification.

Using Salmon
------------

As mentioned above, there are two "modes" of operation for Salmon.  The first,
like Sailfish, requires you to build an index for the transcriptome, but then
subsequently processes reads directly.  The second mode simply requires you to
provide a FASTA file of the transcriptome and a ``.sam`` or ``.bam`` file
containing a set of alignments.

Alignment-based mode
--------------------

Say that you've prepared your alignments using your favorite aligner and the
results are in the file ``aln.bam``, and assume that the sequence of the
transcriptome you want to quantify is in the file ``transcripts.fa``.  You
would run Salmon as follows:

::

    > ./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant

The ``<LIBTYPE>`` parameter is described below and is shared between both modes
of Salmon.  After Salmon has finished running, there will be a directory called
``salmon_quant``, that contains a file called ``quant.sf``.  This contains the
quantification results for the run, and the columns it contains are similar to
those of Sailfish (and self-explanatory where they differ).

For the full set of options that can be passed to Salmon in its alignment-based
mode, and a description of each, run ``salmon quant --help-alignment``.

.. topic:: Multiple alignment files
    
    Currently, salmon assumes that the alignment records for all of the reads 
    corresponding to a single sample appear in the same .bam/.sam file.  This 
    means that if you have multiple .bam/.sam files with aligned reads from the
    *same sample* that you wish to use together for quantification, you must
    first merge the files.  This limitation is temporary and will be removed 
    soon so that multiple different alignment files can be provided for a
    single sample (though they should all be with respect to the same target
    transcripts).


Read-based mode
---------------

If you want to use salmon like sailfish, then you first have to build an salmon
index for your transcriptome.  Again, assume that ``transcripts.fa`` contains
the set of transcripts you wish to quantify.  First, you run the salmon
indexer:

::
    
    > ./bin/salmon index -t transcripts.fa -i transcripts_index

Then, you can quantify any set of reads (say, paired-end reads in files
reads1.fa and reads2.fa) directly against this index using the salmon ``quant``
command as follows:

::

    > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 reads1.fa -2 reads2.fa -o transcripts_quant

You can, of course, pass a number of options to control things such as the
number of threads used or the different cutoffs used for counting reads.
Sust as with the alignment-based mode, after salmon has finished running, there
will be a directory called ``salmon_quant``, that contains a file called
``quant.sf`` containing the quantification results.

What's this ``LIBTYPE``?
------------------------

Salmon, like sailfish, has the user provide a description of the type of
sequencing library from which the reads come, and this contains information
about e.g. the relative orientation of paired end reads.  However, we've
replaced the somewhat esoteric description of the library type with a simple
set of strings; each of which represents a different type of read library. This
new method of specifying the type of read library is being back-ported into
Sailfish and will be available in the next release.

The library type string consists of three parts: the relative orientation of
the reads, the strandedness of the library, and the directionality of the
reads.

The first part of the library string (relative orientation) is only provided if
the library is paired-end. The possible options are:

::

    I = inward
    O = outward
    M = matching

The second part of the read library string specifies whether the protocol is
stranded or unstranded; the options are:

::

    S = stranded
    U = unstranded

If the protocol is unstranded, then we're done.  The final part of the library
string specifies the strand from which the read originates in a strand-specific
protocol — it is only provided if the library is stranded (i.e. if the
library format string is of the form S).  The possible values are:

::

    F = read 1 (or single-end read) comes from the forward strand
    R = read 1 (or single-end read) comes from the reverse strand

An example of some library format strings and their interpretations are:

::

    IU (an unstranded paired-end library where the reads face each other)

::

    SF (a stranded single-end protocol where the reads come from the forward strand)

::

    OSR (a stranded paired-end protocol where the reads face away from each other,
         read1 comes from reverse strand and read2 comes from the forward strand)

Misc
----

Salmon deals with reading from compressed read files in the same way as
sailfish --- by using process substitution.  Say in the read-based salmon
example above, the reads were actually in the files ``reads1.fa.gz`` and
``reads2.fa.gz``, then you'd run the following command to decompress the reads
"on-the-fly":

::

    > ./bin/salmon quant -i transcripts_index -l <LIBTYPE> -1 <(gzcat reads1.fa.gz) -2 <(gzcat reads2.fa.gz) -o transcripts_quant

and the gzipped files will be decompressed via separate processes and the raw
reads will be fed into salmon.

.. note:: Reading through decompressed files multiple times
    Salmon requires a specific number of observations (mapped fragments) to
    be observed before it will report its quantification results.  If it 
    doesn't see enough fragments when reading through the read files the 
    first time, it will read through them again (Don't worry; it's not 
    double counting. The results from the first pass essentially become 
    a "prior" for assigning the proper read counts in subsequent passes).
    However, a named-pipe as created by the process substitution syntax 
    above cannot be read from multiple times.  This means that if your 
    file doesn't have enough mapping fragments you either need to reduce 
    the required number of observations, via the ``-n`` argument, which 
    *may* affect accuracy if it is set too low, or extract the reads to 
    a regular fasta/q file.  We hope to support directly reading from 
    compressed files soon to avoid this necessity.

**Finally**, the purpose of making this beta executable (as well as the Salmon
code) available is for people to use it and provide feedback.  A pre-print and
manuscript are in the works, but the earlier we get feedback, thoughts,
suggestions and ideas, the better!  So, if you have something useful to report
or just some interesting ideas or suggestions, please contact us
(`rob.patro@cs.stonybrook.edu` and/or `carlk@cs.cmu.edu`).  Also, please use
the same e-mail addresses to contact us with any *detailed* bug-reports (though
bug-support for these early beta versions may be slow).
