#ifndef  MINIBATCH_INFO
#define  MINIBATCH_INFO

#include <vector>
#include "LibraryFormat.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include "AlignmentGroup.hpp"

template <typename AlnGroupT>
class MiniBatchInfo {
    public:
        MiniBatchInfo(size_t _batchNum, std::vector<AlnGroupT*>* _alignments, double _logForgettingMass) :
            batchNum(_batchNum), alignments(_alignments), logForgettingMass(_logForgettingMass) {}

        size_t batchNum;
        //AlignmentBatch* alignments;
        //std::vector<size_t>* readGroupBoundaries;
        std::vector<AlnGroupT*>* alignments;
        double logForgettingMass;

        void release(tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue,
                    tbb::concurrent_bounded_queue<AlnGroupT*>& alignmentGroupQueue);
};


template <>
void MiniBatchInfo<AlignmentGroup<ReadPair>>::release(
        tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue,
        tbb::concurrent_bounded_queue<AlignmentGroup<ReadPair>*>& alignmentGroupQueue) {
    size_t ng{0};
    for (auto alnGroup : *alignments) {
        for (auto aln : alnGroup->alignments()) {
            alignmentStructureQueue.push(aln.read1);
            alignmentStructureQueue.push(aln.read2);
            //bam_destroy1(aln.read1); bam_destroy1(aln.read2);
        }
        alnGroup->alignments().clear();
        alignmentGroupQueue.push(alnGroup);
        //delete alnGroup;
        alnGroup = nullptr;
        ++ng;
    }
    delete alignments;
    alignments = nullptr;
}

template <>
void MiniBatchInfo<AlignmentGroup<UnpairedRead>>::release(
        tbb::concurrent_bounded_queue<bam1_t*>& alignmentStructureQueue,
        tbb::concurrent_bounded_queue<AlignmentGroup<UnpairedRead>*>& alignmentGroupQueue) {
    size_t ng{0};
    for (auto alnGroup : *alignments) {
        for (auto aln : alnGroup->alignments()) {
            alignmentStructureQueue.push(aln.read);
            //bam_destroy1(aln.read1); bam_destroy1(aln.read2);
        }
        alnGroup->alignments().clear();
        alignmentGroupQueue.push(alnGroup);
        //delete alnGroup;
        alnGroup = nullptr;
        ++ng;
    }
    delete alignments;
    alignments = nullptr;
}


#endif // MINIBATCH_INFO
