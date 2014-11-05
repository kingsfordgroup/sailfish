#ifndef SALMON_OPTS_HPP
#define SALMON_OPTS_HPP

/**
  * A structure to hold some common options used
  * by Salmon so that we don't have to pass them
  * all around as separate arguments.
  */
struct SalmonOpts {
    // The options below are adopted from the mem_opt_t structure of BWA
    /*
    int maxOccurences; // maximum number of allowable occurences of (S)MEM
    int minSeedOccurences; // try to split a seed into smaller seeds if it occurs
                           // fewer than this many times.
    int minSeedLen; // A seed must be at least this long.
    float splitFactor; // Split a seed if it's longer than splitFactor * minSeedLen.
    int flag; // Used by bwa
    bool maxMEMIntervals; // If true, don't split (S)MEMs into MEMs
    */

    SalmonOpts() : splitSpanningSeeds(false), useFragLenDist(false),
                   useReadCompat(false), maxReadOccs(200) {}
    bool splitSpanningSeeds; // Attempt to split seeds that span multiple transcripts.

    bool useFragLenDist; // Give a fragment assignment a likelihood based on an emperically
                         // observed fragment length distribution.

    bool useReadCompat; // Give a fragment assignment a likelihood based on the compatibility
                        // between the manner in which it mapped and the expected read
                        // librarry format.

    bool useErrorModel; // Learn and apply the error model when computing the likelihood
                        // of a given alignment.

    uint32_t maxReadOccs; // Discard reads  mapping to more than this many places.

    uint32_t maxExpectedReadLen; // Maximum expected length of an observed read.

};

#endif // SALMON_OPTS_HPP
