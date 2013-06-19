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