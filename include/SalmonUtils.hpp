#ifndef __SALMON_UTILS_HPP__
#define __SALMON_UTILS_HPP__

#include "AlignmentLibrary.hpp"
#include "format.h"

#include <boost/filesystem.hpp>

#include <memory>

namespace salmon {
namespace utils {
    template <typename FragT>
    void writeAbundances(AlignmentLibrary<FragT>& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments="");
}
}

#endif // __SALMON_UTILS_HPP__
