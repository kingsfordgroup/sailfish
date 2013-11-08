#ifndef __PARTITION_REFINER_HPP__
#define __PARTITION_REFINER_HPP__

#include <vector>
#include "LookUpTableUtils.hpp"

class PartitionRefiner {
public:
        PartitionRefiner(LUTTools::KmerID numElem);
        void splitWith(std::vector<LUTTools::KmerID> splitter);
        void relabel();
        const std::vector<LUTTools::KmerID>& partitionMembership();

private:
        LUTTools::KmerID numElem_;
        std::vector<LUTTools::KmerID> membership_;
        LUTTools::KmerID maxSetIdx_;
};

#endif // __PARTITION_REFINER_HPP__
