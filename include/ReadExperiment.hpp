#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
}

// Our includes
#include "ClusterForest.hpp"
#include "Transcript.hpp"
#include "ReadLibrary.hpp"

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

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
class ReadExperiment {

    public:

    ReadExperiment(std::vector<ReadLibrary>& readLibraries,
                   //const boost::filesystem::path& transcriptFile,
                   const boost::filesystem::path& indexDirectory) :
        readLibraries_(readLibraries),
        //transcriptFile_(transcriptFile),
        transcripts_(std::vector<Transcript>()),
        batchNum_(0),
        totalAssignedFragments_(0) {
            namespace bfs = boost::filesystem;

            // Make sure the read libraries are valid.
            for (auto& rl : readLibraries_) { rl.checkValid(); }

            // Make sure the transcript file exists.
            /*
            if (!bfs::exists(transcriptFile_)) {
                std::stringstream ss;
                ss << "The provided transcript file: " << transcriptFile_ <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            }
            */

            // ====== Load the transcripts from file
            { // mem-based
                bfs::path indexPath = indexDirectory / "bwaidx";
                if ((idx_ = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS)) == 0) {
                    fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                    fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                    std::exit(1);
                }
            }

            size_t numRecords = idx_->bns->n_seqs;
            std::vector<Transcript> transcripts_tmp;

            fmt::print(stderr, "Index contained {} targets\n", numRecords);
            //transcripts_.resize(numRecords);
            for (auto i : boost::irange(size_t(0), numRecords)) {
                uint32_t id = i;
                char* name = idx_->bns->anns[i].name;
                uint32_t len = idx_->bns->anns[i].len;
                // copy over the length, then we're done.
                transcripts_tmp.emplace_back(id, name, len);
            }

            std::sort(transcripts_tmp.begin(), transcripts_tmp.end(),
                    [](const Transcript& t1, const Transcript& t2) -> bool {
                    return t1.id < t2.id;
                    });

            double alpha = 0.005;
            for (auto& t : transcripts_tmp) {
                transcripts_.emplace_back(t.id, t.RefName.c_str(), t.RefLength, alpha);
            }
            transcripts_tmp.clear();
            // ====== Done loading the transcripts from file

            /*
            FASTAParser fp(transcriptFile.string());

            fmt::print(stderr, "Populating targets from aln = {}, fasta = {} . . .",
                       alignmentFile_, transcriptFile_);
            fp.populateTargets(transcripts_);
            fmt::print(stderr, "done\n");
            */

            // Create the cluster forest for this set of transcripts
            clusters_.reset(new ClusterForest(transcripts_.size(), transcripts_));
        }

    // Free the INDEX
    // ---- Get rid of things we no longer need --------
    //bwa_idx_destroy(idx);

    std::vector<Transcript>& transcripts() { return transcripts_; }

    uint64_t numAssignedFragments() { return numAssignedFragments_; }
    uint64_t numMappedReads() { return numAssignedFragments_; }

    template <typename CallbackT>
    bool processReads(const uint32_t& numThreads, CallbackT& processReadLibrary) {
        bool burnedIn = (totalAssignedFragments_ + numAssignedFragments_ > 5000000);
        for (auto& rl : readLibraries_) {
            processReadLibrary(rl, idx_, transcripts_, clusterForest(),
                               numAssignedFragments_, batchNum_, numThreads, burnedIn);
        }
        return true;
    }

    ClusterForest& clusterForest() { return *clusters_.get(); }

    std::string readFilesAsString() {
        std::stringstream sstr;
        size_t ln{0};
        size_t numReadLibraries{readLibraries_.size()};

        for (auto &rl : readLibraries_) {
            sstr << rl.readFilesAsString();
            if (ln++ < numReadLibraries) { sstr << "; "; }
        }
        return sstr.str();
    }

    bool reset() {
        namespace bfs = boost::filesystem;
        for (auto& rl : readLibraries_) {
            if (!rl.isRegularFile()) { return false; }
        }

        numObservedFragments_ = 0;
        totalAssignedFragments_ += numAssignedFragments_;
        numAssignedFragments_ = 0;
        // batchNum_ = 0; # don't reset batch num right now!
        quantificationPasses_++;
        return true;
    }

    private:
    /**
     * The file from which the alignments will be read.
     * This can be a SAM or BAM file, and can be a regular
     * file or a fifo.
     */
    std::vector<ReadLibrary> readLibraries_;
    /**
     * The file from which the transcripts are read.
     * This is expected to be a FASTA format file.
     */
    //boost::filesystem::path transcriptFile_;
    /**
     * The targets (transcripts) to be quantified.
     */
    std::vector<Transcript> transcripts_;
    /**
     * The index we've built on the set of transcripts.
     */
    bwaidx_t *idx_{nullptr};
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
    std::atomic<uint64_t> numObservedFragments_{0};
    std::atomic<uint64_t> numAssignedFragments_{0};
    std::atomic<uint64_t> batchNum_{0};
    uint64_t totalAssignedFragments_;
    size_t quantificationPasses_;
};

#endif // EXPERIMENT_HPP
