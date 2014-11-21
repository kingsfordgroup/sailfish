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
  * Simple structure holding info about the alignment file.
  */
struct AlignmentFile {
    boost::filesystem::path fileName;
    samFile* fp;
    bam_header_t* header;
};

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
  BAMQueue(std::vector<boost::filesystem::path>& fnames, LibraryFormat& libFmt);
  ~BAMQueue();
  void forceEndParsing();

  bam_header_t* header();
  bam_header_t* safeHeader();

  std::vector<bam_header_t*> headers();

  template <typename FilterT>
  void start(FilterT filt);

  inline bool getAlignmentGroup(AlignmentGroup<FragT*>*& group);

  // Return the number of reads processed so far by the queue
  size_t numObservedReads();
  size_t numMappedReads();

  void reset();

  tbb::concurrent_bounded_queue<FragT*>& getFragmentQueue();
  tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>& getAlignmentGroupQueue();

private:
  size_t popNum{0};
  /** Fill the queue with the appropriate type of alignment
   * depending on the template paramater T
   */
  template <typename FilterT>
  void fillQueue_(FilterT);

  /** Overload of getFrag_ for paired-end reads */
  template <typename FilterT>
  inline bool getFrag_(ReadPair& rpair, FilterT filt);
  /** Overload of getFrag_ for single-end reads */
  template <typename FilterT>
  inline bool getFrag_(UnpairedRead& sread, FilterT filt);

private:
  std::vector<AlignmentFile> files_;
  std::string fname_;
  LibraryFormat libFmt_;

  std::vector<AlignmentFile>::iterator currFile_;
  samFile* fp_ = nullptr;
  bam_header_t* hdr_ = nullptr;

  //htsFile* fp_ = nullptr;
  size_t totalReads_;
  size_t numUnaligned_;
  size_t numMappedReads_;
  tbb::concurrent_bounded_queue<FragT*> fragmentQueue_;
  tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*> alnGroupPool_;
  //tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*> alnGroupQueue_;
  boost::lockfree::spsc_queue<AlignmentGroup<FragT*>*,
                              boost::lockfree::capacity<65535>> alnGroupQueue_;
  volatile bool doneParsing_;
  std::unique_ptr<std::thread> parsingThread_;
  size_t batchNum_;
  std::string readMode_;
};

#include "BAMQueue.tpp"
#endif //__BAMQUEUE_HPP__
