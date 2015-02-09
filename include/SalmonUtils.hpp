#ifndef __SALMON_UTILS_HPP__
#define __SALMON_UTILS_HPP__

extern "C" {
    #include "io_lib/scram.h"
    #include "io_lib/os.h"
}

#include "format.h"

#include "SalmonOpts.hpp"

#include <boost/filesystem.hpp>
#include <memory>

class ReadExperiment;
class LibraryFormat;

namespace salmon {
namespace utils {

    enum class OrphanStatus: uint8_t { LeftOrphan = 0, RightOrphan = 1, Paired = 2 };

    bool headersAreConsistent(SAM_hdr* h1, SAM_hdr* h2);

    bool headersAreConsistent(std::vector<SAM_hdr*>&& headers);

    template <typename AlnLibT>
    void writeAbundances(const SalmonOpts& sopt,
                         AlnLibT& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments="");

    double logAlignFormatProb(const LibraryFormat observed, const LibraryFormat expected, double incompatPrior);

    std::ostream& operator<<(std::ostream& os, OrphanStatus s);
    /**
    *  Given the information about the position and strand from which a paired-end
    *  read originated, return the library format with which it is compatible.
    */
    LibraryFormat hitType(uint32_t end1Start, bool end1Fwd,
                          uint32_t end2Start, bool end2Fwd);
    /**
    *  Given the information about the position and strand from which the
    *  single-end read originated, return the library format with which it
    *  is compatible.
    */
    LibraryFormat hitType(uint32_t readStart, bool isForward);

}
}

#endif // __SALMON_UTILS_HPP__
