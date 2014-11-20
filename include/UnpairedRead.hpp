#ifndef UNPAIRED_READ
#define UNPAIRED_READ

extern "C" {
#ifdef HAVE_CONFIG_H
#undef HAVE_CONFIG_H
#endif

#include "io_lib/scram.h"
#include "io_lib/os.h"

//#define HAVE_CONFIG_H
}

#include "StadenUtils.hpp"
#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"

struct UnpairedRead {
   bam_seq_t* read = nullptr;
   double logProb;

   UnpairedRead() : read(staden::utils::bam_init()), logProb(sailfish::math::LOG_0) {}
   UnpairedRead(bam_seq_t* r, double lp) :
       read(r), logProb(lp) {}

   UnpairedRead(UnpairedRead&& other) {
       logProb = other.logProb;
       std::swap(read, other.read);
   }

   UnpairedRead& operator=(UnpairedRead&& other) {
       logProb = other.logProb;
       std::swap(read, other.read);
       return *this;
   }

   UnpairedRead(const UnpairedRead& other) = default;

   UnpairedRead& operator=(UnpairedRead& other) = default;

   UnpairedRead* clone() {
       return new UnpairedRead(bam_dup(read), logProb);
   }

   ~UnpairedRead() { staden::utils::bam_destroy(read); }

    // return 0 on success, -1 on failure
    int writeToFile(scram_fd* fp) {
        return scram_put_seq(fp, read);
    }

    inline char* getName() {
        return  bam_name(read);
    }

    inline uint32_t getNameLength() {
        return bam_name_len(read);
    }

   inline bool isRight() { return bam_flag(read) & BAM_FREVERSE; }
   inline bool isLeft()  { return !isRight(); }
   inline int32_t left() { return bam_pos(read); }
   inline int32_t right() { return left() + bam_seq_len(read); }
   inline uint32_t fragLen() { return 0; }
   inline ReadType fragType() { return ReadType::SINGLE_END; }
   inline int32_t transcriptID() { return bam_ref(read); }

    inline double logQualProb() {
        int q = bam_map_qual(read);
        //double logP = (q == 255) ? calcQuality(read) : std::log(std::pow(10.0, -q * 0.1));
        double logP = (q == 255) ? sailfish::math::LOG_1 : std::log(std::pow(10.0, -q * 0.1));
        return logP;
    }

};

#endif //UNPAIRED_READ
