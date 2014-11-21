#ifndef READ_PAIR
#define READ_PAIR

extern "C" {
#include "samtools/sam.h"
}

#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"
#include "SalmonUtils.hpp"

#include "format.h"

struct ReadPair {
    bam1_t* read1 = nullptr;
    bam1_t* read2 = nullptr;
    salmon::utils::OrphanStatus orphanStatus;
    double logProb;

    ReadPair():
        read1(bam_init1()), read2(bam_init1()),
        orphanStatus(salmon::utils::OrphanStatus::Paired),
        logProb(sailfish::math::LOG_0) {}


    ReadPair(bam1_t* r1, bam1_t* r2, salmon::utils::OrphanStatus os, double lp) :
        read1(r1), read2(r2), orphanStatus(os), logProb(lp) {}

    ReadPair(ReadPair&& other) {
        orphanStatus = other.orphanStatus;
        logProb = other.logProb;
        std::swap(read1, other.read1);
        std::swap(read2, other.read2);
    }

    ReadPair& operator=(ReadPair&& other) {
        orphanStatus = other.orphanStatus;
        logProb = other.logProb;
        std::swap(read1, other.read1);
        std::swap(read2, other.read2);
        return *this;
    }


   ReadPair(const ReadPair& other) = default;

   ReadPair& operator=(ReadPair& other) = default;

   ReadPair* clone() {
       return new ReadPair(bam_dup1(read1), bam_dup1(read2), orphanStatus, logProb);
   }

    ~ReadPair() {
        bam_destroy1(read1);
        bam_destroy1(read2);
    }

    inline bool isPaired() const { return (orphanStatus == salmon::utils::OrphanStatus::Paired); }
    inline bool isLeftOrphan() const { return (orphanStatus == salmon::utils::OrphanStatus::LeftOrphan); }
    inline bool isRightOrphan() const { return (orphanStatus == salmon::utils::OrphanStatus::RightOrphan); }

    int writeToFile(htsFile* fp) {
        int r1 = bam_write1(fp->fp.bgzf, read1);
        if (r1 > 0 and isPaired()) {
            return bam_write1(fp->fp.bgzf, read2);
        } else {
            return r1;
        }
    }

    inline char* getName() {
        return bam_get_qname(read1);//bam1_qname(read1);
    }

    inline uint32_t getNameLength() {
        return read1->core.l_qname;
    }

    inline uint32_t fragLen() {
        if (!isPaired()) { return 0; }
        auto leftmost1 = read1->core.pos;
        auto leftmost2 = read2->core.pos;

        // The length of the mapped read that is "rightmost" w.r.t. the forward strand.
        auto rightmostLen = (leftmost1 < leftmost2) ? read2->core.l_qseq : read1->core.l_qseq;
        return std::abs(leftmost1 - leftmost2) + rightmostLen;

        //return std::abs(read1->core.isize) + std::abs(read1->core.l_qseq) + std::abs(read2->core.l_qseq);
    }

    inline bool isRight() { return isPaired() ? false : (isRightOrphan() ? true : false) ; }
    inline bool isLeft()  { return isPaired() ? false : (isLeftOrphan() ? true : false); }

    inline int32_t left() {
        if (isPaired()) {
            return std::min(read1->core.pos, read2->core.pos);
        } else {
            return read1->core.pos;
        }
    }

    inline int32_t right() {
        if (isPaired()) {
            return std::max(read1->core.pos + read1->core.l_qseq,
                            read2->core.pos + read2->core.l_qseq);
        } else {
            return read1->core.pos + read1->core.l_qseq;
        }
    }

    inline ReadType fragType() { return ReadType::PAIRED_END; }
    inline int32_t transcriptID() { return read1->core.tid; }

    inline double logQualProb() {
        double q1 = read1->core.qual;
        double q2 = read2->core.qual;
        double logP1 = (q1 == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q1 * 0.1));
        double logP2 = (q2 == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q2 * 0.1));
        return logP1 + logP2;
    }
};

#endif //READ_PAIR
