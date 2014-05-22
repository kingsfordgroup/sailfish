#ifndef UNPAIRED_READ
#define UNPAIRED_READ

extern "C" {
#include "sam.h"
}

#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"

struct UnpairedRead {
   bam1_t* read;
   double logProb;

    inline char* getName() {
        return bam1_qname(read);
    }

   inline bool isRight() { return read->core.flag & BAM_FREVERSE; }
   inline bool isLeft()  { return !isRight(); }
   inline int32_t left() { return read->core.pos; }
   inline int32_t right() { return read->core.pos + read->core.l_qseq; }
   inline uint32_t fragLen() { return 0; }
   inline ReadType fragType() { return ReadType::SINGLE_END; }
   inline int32_t transcriptID() { return read->core.tid; }

    inline double logQualProb() {
        auto q = read->core.qual;
        //double logP = (q == 255) ? calcQuality(read) : std::log(std::pow(10.0, -q * 0.1));
        double logP = (q == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q * 0.1));
        return logP;
    }

};

#endif //UNPAIRED_READ
