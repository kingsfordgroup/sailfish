#ifndef ALIGNMENT_LIBRARY_HPP
#define ALIGNMENT_LIBRARY_HPP

extern "C" {
#include "htslib/sam.h"
#include "samtools/samtools.h"
}

// Our includes
#include "Transcript.hpp"
#include "BAMQueue.hpp"
#include "LibraryFormat.hpp"
#include "AlignmentGroup.hpp"
#include "FASTAParser.hpp"

// Boost includes
#include <boost/filesystem.hpp>

// Standard includes
#include <vector>
#include <memory>

/**
  *  This class represents an library of alignments used to quantify
  *  a set of target transcripts.  The AlignmentLibrary contains info
  *  about both the alignment file and the target sequence (transcripts).
  *  It is used to group them together and track information about them
  *  during the quantification procedure.
  */
template <typename FragT>
class AlignmentLibrary {

    public:

    AlignmentLibrary(const boost::filesystem::path& alnFile,
                     const boost::filesystem::path& transcriptFile,
                     LibraryFormat libFmt) :
        alignmentFile_(alnFile),
        transcriptFile_(transcriptFile),
        libFmt_(libFmt),
        transcripts_(std::vector<Transcript>()),
        quantificationPasses_(0) {
            namespace bfs = boost::filesystem;

            // Make sure the alignment file exists.
            if (!bfs::exists(alignmentFile_)) {
                std::stringstream ss;
                ss << "The provided alignment file: " << alignmentFile_ <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            }

            // Make sure the transcript file exists.
            if (!bfs::exists(transcriptFile_)) {
                std::stringstream ss;
                ss << "The provided transcript file: " << transcriptFile_ <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            }

            // The alignment file existed, so create the alignment queue
            std::string fname = alignmentFile_.string();
            bq = std::unique_ptr<BAMQueue<FragT>>(new BAMQueue<FragT>(fname, libFmt_));
            bam_header_t* header = bq->header();

            // The transcript file existed, so load up the transcripts
            double alpha = 0.005;
            for (size_t i = 0; i < header->n_targets; ++i) {
                transcripts_.emplace_back(i, header->target_name[i], header->target_len[i], alpha);
            }

            FASTAParser fp(transcriptFile.string());

            fmt::print(stderr, "Populating targets from aln = {}, fasta = {} . . .",
                    alignmentFile_, transcriptFile_);
            fp.populateTargets(transcripts_);
            fmt::print(stderr, "done\n");

            // Start parsing the alignments
            bq->start();
        }

    std::vector<Transcript>& transcripts() { return transcripts_; }

    inline bool getAlignmentGroup(AlignmentGroup<FragT>*& ag) { return bq->getAlignmentGroup(ag); }

    inline tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue() {
        return bq->getAlignmentStructureQueue();
    }

    inline tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*>& alignmentGroupQueue() {
        return bq->getAlignmentGroupQueue();
    }

    inline size_t numMappedReads() { return bq->numMappedReads(); }
    const boost::filesystem::path& alignmentFile() { return alignmentFile_; }

    bool reset() {
        namespace bfs = boost::filesystem;
        if (!bfs::is_regular_file(alignmentFile_)) {
            return false;
        }
        std::string fname = alignmentFile_.string();
        bq->reset();
        bq->start();
        quantificationPasses_++;
        return true;
    }

    private:
    boost::filesystem::path alignmentFile_;
    boost::filesystem::path transcriptFile_;
    LibraryFormat libFmt_;
    std::vector<Transcript> transcripts_;
    std::unique_ptr<BAMQueue<FragT>> bq;
    size_t quantificationPasses_;
};

#endif // ALIGNMENT_LIBRARY_HPP
