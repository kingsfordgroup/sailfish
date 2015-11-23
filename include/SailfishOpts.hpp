#ifndef __SAILFISH_OPTS_HPP__
#define __SAILFISH_OPTS_HPP__

#include "spdlog/spdlog.h"
#include <boost/filesystem.hpp>
#include <memory>

struct SailfishOpts {
    uint32_t numThreads; // number of threads to use
    bool allowOrphans;
    bool noEffectiveLengthCorrection;
    bool noFragLengthDist;
    bool useVBOpt{false};
    bool useGSOpt{false};
    bool useUnsmoothedFLD{false};
    bool ignoreLibCompat{false};
    bool enforceLibCompat{false};
    bool allowDovetail{false};
    uint32_t maxFragLen;
    uint32_t numGibbsSamples;
    uint32_t numBootstraps;
    uint32_t maxReadOccs;
    size_t fragLenDistMax;
    size_t fragLenDistPriorMean;
    size_t fragLenDistPriorSD;
    std::shared_ptr<spdlog::logger> jointLog{nullptr};
    std::shared_ptr<spdlog::logger> fileLog{nullptr};
    boost::filesystem::path indexDirectory;
    boost::filesystem::path outputDirectory;
};

#endif // __SAILFISH_OPTS_HPP__
