#ifndef COLLAPSED_EM_OPTIMIZER_HPP
#define COLLAPSED_EM_OPTIMIZER_HPP

#include <memory>
#include <unordered_map>
#include <functional>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "ReadExperiment.hpp"
#include "SailfishOpts.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

class BootstrapWriter;

class CollapsedEMOptimizer {
    public:
        using VecType = std::vector<tbb::atomic<double>>;
        using SerialVecType = std::vector<double>;
        CollapsedEMOptimizer();

        bool optimize(ReadExperiment& readExp,
                      SailfishOpts& sopt,
                      double tolerance = 0.01,
                      uint32_t maxIter = 1000);

        bool gatherBootstraps(
                ReadExperiment& readExp,
                SailfishOpts& sopt,
	        std::function<bool(const std::vector<double>&)>& writeBootstrap,
                double relDiffTolerance,
                uint32_t maxIter);
};

#endif // COLLAPSED_EM_OPTIMIZER_HPP

