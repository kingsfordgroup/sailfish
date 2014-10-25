#ifndef READ_PAIR
#define READ_PAIR

extern "C" {
#include "samtools/sam.h"
}

#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"

#include "format.h"

struct ReadPair {
    bam1_t* read1 = nullptr;
    bam1_t* read2 = nullptr;
    double logProb;

    ReadPair():
        read1(bam_init1()), read2(bam_init1()), logProb(sailfish::math::LOG_0) {}


    ReadPair(bam1_t* r1, bam1_t* r2, double lp) :
        read1(r1), read2(r2), logProb(lp) {}

    ReadPair(ReadPair&& other) {
        logProb = other.logProb;
        std::swap(read1, other.read1);
        std::swap(read2, other.read2);
    }

    ReadPair& operator=(ReadPair&& other) {
        logProb = other.logProb;
        std::swap(read1, other.read1);
        std::swap(read2, other.read2);
        return *this;
    }

   ReadPair(const ReadPair& other) = default;

   ReadPair& operator=(ReadPair& other) = default;


    ~ReadPair() {
        bam_destroy1(read1);
        bam_destroy1(read2);
    }


    inline char* getName() {
        return bam_get_qname(read1);//bam1_qname(read1);
    }

    inline uint32_t getNameLength() {
        return read1->core.l_qname;
    }

    inline uint32_t fragLen() {
        //return std::abs(read1->core.pos - read2->core.pos) + read2->core.l_qseq;
        return std::abs(read1->core.isize) + std::abs(read1->core.l_qseq) + std::abs(read2->core.l_qseq);
    }

    inline bool isRight() { return false; }
    inline bool isLeft()  { return false; }

    inline int32_t left() { return std::min(read1->core.pos, read2->core.pos); }
    inline int32_t right() {
        return std::max(read1->core.pos + read1->core.l_qseq,
                        read2->core.pos + read2->core.l_qseq);
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
