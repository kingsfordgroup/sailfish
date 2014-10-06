#ifndef __SALMON_UTILS_HPP__
#define __SALMON_UTILS_HPP__

#include "AlignmentLibrary.hpp"
#include "format.h"

#include <boost/filesystem.hpp>

#include <memory>

class ReadExperiment;
class LibraryFormat;

namespace salmon {
namespace utils {
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
