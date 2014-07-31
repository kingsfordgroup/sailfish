#ifndef ALIGNMENT_GROUP
#define ALIGNMENT_GROUP

extern "C" {
#include "htslib/sam.h"
#include "samtools/samtools.h"
}

#include <vector>
#include "SailfishMath.hpp"
#include "ReadPair.hpp"

template <typename FragT>
class AlignmentGroup {
    public:
        AlignmentGroup(){ alignments_.reserve(10); }
        std::vector<FragT>& alignments() { return alignments_; }
        void addAlignment(const FragT& p) { alignments_.push_back(p);}
        void addAlignment(bam1_t* r) {
            alignments_.push_back({r, sailfish::math::LOG_0});
        }
        void addAlignment(bam1_t* r1, bam1_t* r2) {
            alignments_.push_back({r1, r2, sailfish::math::LOG_0});
        }
        size_t size() { return alignments_.size(); }
    private:
        std::vector<FragT> alignments_;
};
#endif // ALIGNMENT_GROUP
