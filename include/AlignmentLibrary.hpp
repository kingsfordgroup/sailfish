#ifndef ALIGNMENT_LIBRARY_HPP
#define ALIGNMENT_LIBRARY_HPP

// samtools / htslib includes
extern "C" {
#include "htslib/sam.h"
#include "samtools/samtools.h"
}

// Our includes
#include "ClusterForest.hpp"
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
  *  This class represents a library of alignments used to quantify
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

            // Create the cluster forest for this set of transcripts
            clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));

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

    inline BAMQueue<FragT>& getAlignmentGroupQueue() { return *bq.get(); }

    inline size_t numMappedReads() { return bq->numMappedReads(); }
    const boost::filesystem::path& alignmentFile() { return alignmentFile_; }
    ClusterForest& clusterForest() { return *clusters_.get(); }

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
    /**
     * The file from which the alignments will be read.
     * This can be a SAM or BAM file, and can be a regular
     * file or a fifo.
     */
    boost::filesystem::path alignmentFile_;
    /**
     * The file from which the transcripts are read.
     * This is expected to be a FASTA format file.
     */
    boost::filesystem::path transcriptFile_;
    /**
     * Describes the expected format of the sequencing
     * fragment library.
     */
    LibraryFormat libFmt_;
    /**
     * The targets (transcripts) to be quantified.
     */
    std::vector<Transcript> transcripts_;
    /**
     * A pointer to the queue from which the fragments
     * will be read.
     */
    std::unique_ptr<BAMQueue<FragT>> bq;
    /**
     * The cluster forest maintains the dynamic relationship
     * defined by transcripts and reads --- if two transcripts
     * share an ambiguously mapped read, then they are placed
     * in the same cluster.
     */
    std::unique_ptr<ClusterForest> clusters_;
    /** Keeps track of the number of passes that have been
     *  made through the alignment file.
     */
    size_t quantificationPasses_;
};

#endif // ALIGNMENT_LIBRARY_HPP
