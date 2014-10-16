.. _FragLibType:

Fragment Library Types
======================

There are numerous library preparation protocols for RNA-seq that result in
sequencing reads with different characteristics.  For example, reads can be
single end (only one side of a fragment is recorded as a read) or paired-end
(reads are generated from both ends of a fragment).  Further, the sequencing
reads themselves may be unstraned or strand-specific.  Finally, paired-end
protocols will have a specified relative orientation.  To characterize the
various different typs of sequencing libraries, we've created a miniature
"language" that allows for the succinct description of the many different types
of possible fragment libraries.  For paired-end reads, the possible
orientations, along with a graphical description of what they mean, are
illustrated below:

.. image:: ReadLibraryIllustration.png

So, for example, if you wanted to specify a fragment library of strand-specific
paired-end reads, oriented toward each other, where read 1 comes from the
forward strand and read 2 comes from the reverse strand, you would specify ``-l
ISF`` on the command line.  This designates that the library being processed has
the type "ISF" meaning, **I**\ nward (the relative orientation), **S**\ tranted
(the protocol is strand-specific), **F**\ orward (read 1 comes from the forward
strand).

The single end library strings are a bit simpler than their pair-end counter
parts, since there is no relative orientation of which to speak.  Thus, the
only possible library format types for single-end reads are ``U`` (for
unstranded), ``SF`` (for strand-specific reads coming from the forward strand)
and ``SR`` (for strand-specific reads coming from the reverse strand).
