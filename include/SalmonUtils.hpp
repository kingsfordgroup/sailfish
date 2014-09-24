#ifndef __SALMON_UTILS_HPP__
#define __SALMON_UTILS_HPP__

#include "AlignmentLibrary.hpp"
#include "format.h"

#include <boost/filesystem.hpp>

#include <memory>

class ReadExperiment;

namespace salmon {
namespace utils {
    template <typename AlnLibT>
    void writeAbundances(AlnLibT& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments="");
}
}

#endif // __SALMON_UTILS_HPP__
