#ifndef EXPERIMENT_HPP
#define EXPERIMENT_HPP

// Our includes
#include "Transcript.hpp"
#include "ReadLibrary.hpp"
//#include "FragmentLengthDistribution.hpp"
//#include "SequenceBiasModel.hpp"
#include "SailfishOpts.hpp"
#include "SailfishIndex.hpp"
#include "EquivalenceClassBuilder.hpp"

// Logger includes
#include "spdlog/spdlog.h"

// Boost includes
#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

// Standard includes
#include <vector>
#include <memory>
#include <fstream>

//S_AYUSH_CODE
#include "ReadKmerDist.hpp"
//T_AYUSH_CODE

/**
  *  This class represents a library of reads used to quantify
  *  a set of target transcripts.
  */
class ReadExperiment {

    public:

    ReadExperiment(std::vector<ReadLibrary>& readLibraries,
                   const boost::filesystem::path& indexDirectory,
		           SailfishOpts& sopt) :
        readLibraries_(readLibraries),
        transcripts_(std::vector<Transcript>()),
        //fragStartDists_(5),
        //seqBiasModel_(1.0),
    	eqBuilder_(sopt.jointLog) {
            namespace bfs = boost::filesystem;

            // Make sure the read libraries are valid.
            for (auto& rl : readLibraries_) { rl.checkValid(); }

            size_t maxFragLen = sopt.fragLenDistMax;
            size_t meanFragLen = sopt.fragLenDistPriorMean;
            size_t fragLenStd = sopt.fragLenDistPriorSD;

            sfIndex_.reset(new SailfishIndex(sopt.jointLog));
            sfIndex_->load(indexDirectory);
            loadTranscriptsFromQuasi();
        }

    EquivalenceClassBuilder& equivalenceClassBuilder() {
        return eqBuilder_;
    }

    std::vector<Transcript>& transcripts() { return transcripts_; }

    const std::vector<Transcript>& transcripts() const { return transcripts_; }

    uint64_t numFragHits() { return numFragHits_; }
    std::atomic<uint64_t>& numFragHitsAtomic() { return numFragHits_; }

    uint64_t numMappedFragments() { return numMappedFragments_; }

    uint64_t upperBoundHits() { return upperBoundHits_; }
    std::atomic<uint64_t>& upperBoundHitsAtomic() { return upperBoundHits_; }
    void setUpperBoundHits(uint64_t ubh) { upperBoundHits_.store(ubh); }

    uint64_t numObservedFragments() const { return numObservedFragments_; }
    std::atomic<uint64_t>& numObservedFragmentsAtomic() { return numObservedFragments_; }

    std::atomic<uint64_t>& numMappedFragmentsAtomic() { return numMappedFragments_; }

    void setNumObservedFragments(uint64_t numObserved) { numObservedFragments_ = numObserved; }

    double mappingRate() {
        return static_cast<double>(numMappedFragments_) / numObservedFragments_;
    }

    SailfishIndex* getIndex() { return sfIndex_.get(); }

    template <typename IndexT>
    void loadTranscriptsFromQuasiIndex(IndexT* idx_) {
        size_t numRecords = idx_->txpNames.size();

        fmt::print(stderr, "Index contained {} targets\n", numRecords);
        double alpha = 0.005;
        for (auto i : boost::irange(size_t(0), numRecords)) {
            uint32_t id = i;
            const char* name = idx_->txpNames[i].c_str();
            uint32_t len = idx_->txpLens[i];
            // copy over the length, then we're done.
            transcripts_.emplace_back(id, name, len);
            auto& txp = transcripts_.back();
            // The transcript sequence
            auto txpSeq = idx_->seq.substr(idx_->txpOffsets[i], len);
            txp.Sequence = txpSeq;
        }
        // ====== Done loading the transcripts from file
    }

    void loadTranscriptsFromQuasi() {
        if (sfIndex_->is64BitQuasi()) {
            loadTranscriptsFromQuasiIndex<RapMapSAIndex<int64_t>>(
                    sfIndex_->quasiIndex64());
        } else {
            loadTranscriptsFromQuasiIndex<RapMapSAIndex<int32_t>>(
                    sfIndex_->quasiIndex32());
        }
	}

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

    double effectiveMappingRate() {
        return effectiveMappingRate_;
    }

    void setEffetiveMappingRate(double emr) {
        effectiveMappingRate_ = emr;
    }
    /*
    SequenceBiasModel& sequenceBiasModel() {
        return seqBiasModel_;
    }
    */

    std::vector<ReadLibrary>& readLibraries() { return readLibraries_; }
    //FragmentLengthDistribution* fragmentLengthDistribution() { return fragLengthDist_.get(); }

    //S_AYUSH_CODE
    ReadKmerDist<6, std::atomic<uint32_t>>& readBias() { return readBias_; }
    //T_AYUSH_CODE

    private:
    /**
     * The file from which the alignments will be read.
     * This can be a SAM or BAM file, and can be a regular
     * file or a fifo.
     */
    std::vector<ReadLibrary> readLibraries_;
    /**
     * The targets (transcripts) to be quantified.
     */
    std::vector<Transcript> transcripts_;
    /**
     * The index we've built on the set of transcripts.
     */
    std::unique_ptr<SailfishIndex> sfIndex_{nullptr};

    //SequenceBiasModel seqBiasModel_;

    /** Keeps track of the number of passes that have been
     *  made through the alignment file.
     */
    std::atomic<uint64_t> numObservedFragments_{0};
    std::atomic<uint64_t> numMappedFragments_{0};
    std::atomic<uint64_t> numFragHits_{0};
    std::atomic<uint64_t> upperBoundHits_{0};
    double effectiveMappingRate_{0.0};
    //std::unique_ptr<FragmentLengthDistribution> fragLengthDist_;
    EquivalenceClassBuilder eqBuilder_;

    //S_AYUSH_CODE
    // Since multiple threads can touch this dist, we
    // need atomic counters.
    ReadKmerDist<6, std::atomic<uint32_t>> readBias_;
    //T_AYUSH_CODE
};

#endif // EXPERIMENT_HPP
