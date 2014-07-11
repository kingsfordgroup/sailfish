#ifndef __BAMQUEUE_HPP__
#define __BAMQUEUE_HPP__

extern "C" {
#include "htslib/sam.h"
#include "samtools/samtools.h"
}

#include <boost/lockfree/spsc_queue.hpp>
#include <boost/lockfree/queue.hpp>
#include <tbb/atomic.h>
#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <memory>
#include <exception>
#include <thread>
#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <tbb/concurrent_queue.h>
#include "AlignmentGroup.hpp"
#include "LibraryFormat.hpp"
#include "SailfishMath.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"

/**
  * A queue from which to draw BAM alignments.  The queue is thread-safe, and
  * can be written to and read from multiple threads.
  *
  * This class is templated on LibT --- the type of read library from which
  * the provided alignments are being generated.
  */
template <typename FragT>
class BAMQueue {
public:
  BAMQueue(const std::string& fname, LibraryFormat& libFmt);
  ~BAMQueue();
  void forceEndParsing();

  bam_header_t* header();

  void start();

  inline bool getAlignmentGroup(AlignmentGroup<FragT>*& group);


  tbb::concurrent_bounded_queue<bam1_t*>& getAlignmentStructureQueue();
  tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*>& getAlignmentGroupQueue();

private:
  size_t popNum{0};
  /** Fill the queue with the appropriate type of alignment
   * depending on the template paramater T
   */
  void fillQueue_();

  /** Overload of getFrag_ for paired-end reads */
  inline bool getFrag_(ReadPair& rpair);
  /** Overload of getFrag_ for single-end reads */
  inline bool getFrag_(UnpairedRead& sread);

private:
  std::string fname_;
  LibraryFormat libFmt_;
  samFile* fp_ = nullptr;
  bam_hdr_t* hdr_ = nullptr;
  //htsFile* fp_ = nullptr;
  size_t totalReads_;
  size_t numUnaligned_;
  tbb::concurrent_bounded_queue<bam1_t*> alnStructQueue_;
  tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*> alnGroupPool_;
  boost::lockfree::spsc_queue<AlignmentGroup<FragT>*,
                              boost::lockfree::capacity<65535>> alnGroupQueue_;
  bool doneParsing_;
  std::thread* parsingThread_ = nullptr;
  size_t batchNum_;
  double logForgettingMass_;
  double forgettingFactor_;
};

#include "BAMQueue.tpp"
#endif //__BAMQUEUE_HPP__
