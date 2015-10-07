Sailfish
================

Sailfish is a tool for transcript quantification from RNA-seq data.  It
requires a set of target transcripts (either from a reference or *de-novo*
assembly) to quantify.  All you need to run sailfish is a fasta file containing
your reference transcripts and a (set of) fasta/fastq file(s) containing your
reads.  Sailfish runs in two phases; indexing and quantification.  The indexing
step is independent of the reads, and only needs to be run once for a particular
set of reference transcripts and choice of k (the k-mer size). The
quantification step, obviously, is specific to the set of RNA-seq reads and is
thus run more frequently.

Indexing
--------

To generate the sailfish index for your reference set of transcripts, you
should run the following command:

::

    > sailfish index -t <ref_transcripts> -o <out_dir> -k <kmer_len>


This will build a sailfish index using k-mers of length ``<kmer_len>`` for the
reference transcripts  provided in the file ``<ref_transcripts>`` and place the
index under the directory ``<out_dir>``.  There  are additional options that can
be passed to the sailfish indexer (e.g. the number of threads to use).  These
can be seen by executing the command ``sailfish index -h``.

Note that, as of v0.7.0, the meaning of the ``-k`` parameter has changed slightly.
Rather than the k-mer size on which Sailfish will quantify abundances, it becomes
the minimum match size that will be considered in the `quasi-mapping <http://github.com/COMBINE-lab/RapMap>`_
procedure during quantification.  For sufficiently long (e.g. 75bp or greater)
reads, the default should be acceptable.  If you have substantially shorter
reads, you may want to consider a smaller ``-k``.

.. note:: values of k

  The ``k`` value used to build the Sailfish index must be an odd number.  Using an
  even value for ``k`` will raise an error and the full index will not be built.

Quantification
--------------

Now that you have generated the sailfish index (say that it's the directory
``<index_dir>`` --- this corresponds to the <out_dir> argument provided in the
previous step), you can quantify the transcript expression for a given set of
reads.  To perform the quantification, you run a command like the following:

::

    > sailfish quant -i <index_dir> -l "<libtype>" {-r <unmated> | -1 <mates1> -2 <mates2>} -o <quant_dir>

Where ``<index_dir>`` is, as described above, the location of the sailfish
index, ``<libtype>`` is a string describing the format of the fragment (read)
library (see :ref:`FragLibType`), ``<unmated>`` is a list of files
containing unmated reads, ``<mates{1,2}>`` are lists of files containg,
respectively, the first and second mates of paired-end reads. Finally,
``<quant_dir>`` is the directory where the output should be written. Just like the
indexing step, additional options are available, and can be viewed by running
``sailfish quant -h``.

When the quantification step is finished, the directory ``<quant_dir>`` will
contain a file named "quant.sf" (and, if bias correction is enabled, an
additional file names "quant_bias_corrected.sf").  This file contains the
result of the Sailfish quantification step.  This file contains a number of
columns (which are listed in the last of the header lines beginning with '#').
Specifically, the columns are (1) Transcript ID, (2) Transcript Length, (3)
Transcripts per Million (TPM) and (6) Estimated number of reads (an estimate
of the number of reads drawn from this transcript given the transcript's
relative abundance and length). The first two columns are self-explanatory,
the next four are measures of transcript abundance and the final is a commonly
used input for differential expression tools.  The Transcripts per Million
quantification number is computed as described in [1]_, and is meant as an
estimate of the number of transcripts, per million observed transcripts,
originating from each isoform.  Its benefit over the F/RPKM measure is that it
is independent of the mean expressed transcript length (i.e. if the mean
expressed transcript length varies between samples, for example, this alone can
affect differential analysis based on the K/RPKM.).

Description of important options
--------------------------------

Sailfish exposes a number of useful optional command-line parameters to the user.
The particularly important ones are explained here, but you can always run
``sailfish quant -h`` to see them all.

""""""""""""""""""""""""""
``-p`` / ``--numThreads``
""""""""""""""""""""""""""

The number of threads that will be used for quasi-mapping, quantification, and
bootstrapping / posterior sampling (if enabled).  Sailfish is designed to work
well with many threads, so, if you have a sufficient number of processors, larger
values here can speed up the run substantially.


""""""""""""""
``--useVBOpt``
""""""""""""""

Use the variational Bayesian EM algorithm rather than the "standard" EM algorithm
to optimize abundance estimates.  The details of the VBEM algorithm can be found
in [2]_, and the details of the variant over fragment equivalence classes that
we use can be found in [3]_.  While both the standard EM and the VBEM produce
accurate abundance estimates, those produced by the VBEM seem, generally, to be
a bit more accurate.  Further, the VBEM tends to converge after fewer iterations,
so it may result in a shorter runtime; especially if you are computing many
bootstrap samples. 

"""""""""""""""""""
``--numBootstraps``
"""""""""""""""""""

Sailfish has the ability to optionally compute bootstrapped abundance estimates.
This is done by resampling (with replacement) from the counts assigned to
the fragment equivalence classes, and then re-running the optimization procedure,
either the EM or VBEM, for each such sample.  The values of these different
bootstraps allows us to assess technical variance in the main abundance estimates
we produce.  Such estimates can be useful for downstream (e.g. differential
expression) tools that can make use of such uncertainty estimates.  This option
takes a positive integer that dictates the number of bootstrap samples to compute.
The more samples computed, the better the estimates of varaiance, but the
more computation (and time) required.

"""""""""""""""""""""
``--numGibbsSamples``
"""""""""""""""""""""

Just as with the bootstrap procedure above, this option produces samples that allow
us to estimate the variance in abundance estimates.  However, in this case the
samples are generated using posterior Gibbs sampling over the fragment equivalence
classes rather than bootstrapping.  We are currently analyzing these different approaches
to assess the potential trade-offs in time / accuracy.  The ``--numBootstraps`` and
``--numGibbsSamples`` options are mutually exclusive (i.e. in a given run, you must
set at most one of these options to a positive integer.)

References
----------

.. [1] Li, Bo, et al. "RNA-Seq gene expression estimation with read mapping uncertainty."
    Bioinformatics 26.4 (2010): 493-500.
.. [2] Nariai, Naoki, et al. "TIGAR: transcript isoform abundance estimation method with gapped alignment of RNA-Seq data by variational Bayesian inference."
    Bioinformatics (2013): btt381.
.. [3] Rob Patro, Geet Duggal & Carl Kingsford "Accurate, fast, and model-aware transcript expression quantification with Salmon"
    bioRxiv doi: http://dx.doi.org/10.1101/021592

.. _CMake : http://www.cmake.org
.. _Boost: http://www.boost.org
