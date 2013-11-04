Sailfish Changelog
===================

v0.1.0
------
* Initial release

v0.1.1
------

__Date__: June 19, 2013

* Speed improvements in "quant" phase -- Substanially
  improved the speed of direction-aware kmer counting.
  In most cases, it should now be almost as fast as counting
  kmers with the "canonical" flag set.  Also improved
  the speed of the EM proceedure by consildating non-dependent
  computations into the same loops so that fewer passes
  over the data are required.

* Fixed bug in computation of the log-likelihood function.
  The SQUAREM procedure now acts upon (and prints out) the
  expected per-transcript log-likelihood rather than the 
  log-likelihood for all transcripts.  The value is more 
  reasonable and works better for the SQUAREM step-size 
  selection procedure.

* Build improvements
    * Removed build dependency on libprofiler 
      (_thanks to Chirs Hill for mentioning this dependency_).
    * Fixed std::this_thread dependence that would lead to
      compilation problems on versions of GCC built without
      "--enable-libstdcxx-time"
      (_thanks to Chirs Hill for the suggested fix_).

v0.5.0
-------

* Initial public release

v0.6.0
-------

* Moved computation of k-mer equivalence classes to the index-building
  phase.  This substantially reduces the memory usage during estimation
  as well as the size of several of the stored indexes.  The algorithm
  used to compute the equivalence classes was also changed from a 
  parallel-hashing based algorithm to a divisive partition refinement 
  algorithm.  This latter algorithm is more suitable to the per-transcript
  processing that happens during the indexing phase.

* Implemented reading from named pipes and input redirection.  Sequencing
  reads can now be streamed in from a named pipe (e.g. using process substitution
  syntax.)

* Implemented direction-aware k-mer counting.  If the directionality (sense / anti-sense)
  for a set of reads is known (e.g. as the result of a direction-aware protocol), it can
  now be specified on the command line.  Thus, there are conceptually 3 "classes" of reads;
  forward/sense, reverse/anti-sense and undirected.

* The estimated number of reads originating from each transcript is now written to the output
  file.  This may be useful for differential-expression tools which are based on read counts.

* Fixed oversight in bias-correction phase where only RPKM estimates (and not e.g. TPMs) 
  were corrected.  Now all different estimates are corrected during the bias-correction
  phase.
