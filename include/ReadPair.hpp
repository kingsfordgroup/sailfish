#ifndef READ_PAIR
#define READ_PAIR

#include "StadenUtils.hpp"
#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"

#include "format.h"

struct ReadPair {
    bam_seq_t* read1 = nullptr;
    bam_seq_t* read2 = nullptr;
    double logProb;

    ReadPair():
        read1(staden::utils::bam_init()), read2(staden::utils::bam_init()), logProb(sailfish::math::LOG_0) {}

    ReadPair(bam_seq_t* r1, bam_seq_t* r2, double lp) :
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

   ReadPair* clone() {
       return new ReadPair(bam_dup(read1), bam_dup(read2), logProb);
   }

    ~ReadPair() {
        staden::utils::bam_destroy(read1);
        staden::utils::bam_destroy(read2);
    }

    /**
      * returns 0 on success, -1 on failure.
      */
    int writeToFile(scram_fd* fp) {
        int r1 = scram_put_seq(fp, read1);
        if (r1 == 0) {
            return scram_put_seq(fp, read2);
        } else {
            return r1;
        }
    }

    inline char* getName() {
        return bam_name(read1);//bam1_qname(read1);
    }

    inline uint32_t getNameLength() {
        return bam_name_len(read1);
    }

    inline uint32_t fragLen() {
        auto leftmost1 = bam_pos(read1);
        auto leftmost2 = bam_pos(read2);

        // The length of the mapped read that is "rightmost" w.r.t. the forward strand.
        auto rightmostLen = (leftmost1 < leftmost2) ? bam_seq_len(read2) : bam_seq_len(read1);
        return std::abs(leftmost1 - leftmost2) + rightmostLen;

        //return std::abs(read1->core.isize) + std::abs(read1->core.l_qseq) + std::abs(read2->core.l_qseq);
    }

    inline bool isRight() { return false; }
    inline bool isLeft()  { return false; }

    inline int32_t left() { return std::min(bam_pos(read1), bam_pos(read2)); }
    inline int32_t right() {
        return std::max(bam_pos(read1) + bam_seq_len(read1),
                        bam_pos(read2) + bam_seq_len(read2));
    }

    inline ReadType fragType() { return ReadType::PAIRED_END; }
    inline int32_t transcriptID() { return bam_ref(read1); }

    inline double logQualProb() {
        double q1 = bam_map_qual(read1);
        double q2 = bam_map_qual(read2);
        double logP1 = (q1 == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q1 * 0.1));
        double logP2 = (q2 == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q2 * 0.1));
        return logP1 + logP2;
    }
};

#endif //READ_PAIR
