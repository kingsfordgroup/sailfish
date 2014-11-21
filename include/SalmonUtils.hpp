#ifndef __SALMON_UTILS_HPP__
#define __SALMON_UTILS_HPP__

extern "C" {
    #include "samtools/bam.h"
}

#include "format.h"

#include <boost/filesystem.hpp>

#include <memory>

class ReadExperiment;
class LibraryFormat;

namespace salmon {
namespace utils {

    enum class OrphanStatus: uint8_t { LeftOrphan = 0, RightOrphan = 1, Paired = 2 };

    bool headersAreConsistent(bam_header_t* h1, bam_header_t* h2);

    bool headersAreConsistent(std::vector<bam_header_t*>&& headers);

    template <typename AlnLibT>
    void writeAbundances(AlnLibT& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments="");
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
