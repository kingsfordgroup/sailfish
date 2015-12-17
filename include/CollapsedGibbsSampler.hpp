#ifndef COLLAPSED_GIBBS_SAMPLER_HPP
#define COLLAPSED_GIBBS_SAMPLER_HPP

#include <functional>
#include <unordered_map>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#define __IS_SAILFISH__

#ifdef __IS_SAILFISH__
#include "SailfishOpts.hpp"
#else
#include "SalmonOpts.hpp"
#endif

#include "cuckoohash_map.hh"

class BootstrapWriter;

class CollapsedGibbsSampler {
    public:
        using VecType = std::vector<double>;
        CollapsedGibbsSampler();

        template <typename ExpT>
        bool sample(ExpT& readExp,
                    SailfishOpts& sopt,
                    std::function<bool(const std::vector<int>&)>& writeSample, 
                    uint32_t numSamples = 500);

};

#endif // COLLAPSED_EM_OPTIMIZER_HPP

