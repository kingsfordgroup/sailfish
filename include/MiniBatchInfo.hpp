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

        template <typename FragT>
        void release(tbb::concurrent_bounded_queue<FragT*>& fragmentQueue,
                    tbb::concurrent_bounded_queue<AlnGroupT*>& alignmentGroupQueue){
            size_t ng{0};
            for (auto& alnGroup : *alignments) {
                for (auto aln : alnGroup->alignments()) {
                    fragmentQueue.push(aln);
                    aln = nullptr;
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



};

/*
template <>
void MiniBatchInfo<AlignmentGroup<ReadPair>>::release(
        tbb::concurrent_bounded_queue<ReadPair*>& alignmentStructureQueue,
        tbb::concurrent_bounded_queue<AlignmentGroup<ReadPair>*>& alignmentGroupQueue) {
    size_t ng{0};
    for (auto& alnGroup : *alignments) {
        for (auto& aln : alnGroup->alignments()) {

            alignmentStructureQueue.push(aln.read1);
            alignmentStructureQueue.push(aln.read2);
            aln.read1 = nullptr;
            aln.read2 = nullptr;
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
*/

#endif // MINIBATCH_INFO
