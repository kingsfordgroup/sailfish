#ifndef __PARTITION_REFINER_HPP__
#define __PARTITION_REFINER_HPP__

#include <vector>

class PartitionRefiner {
public:
        PartitionRefiner(std::size_t numElem);
        void splitWith(std::vector<std::size_t> splitter);
        void relabel();
        const std::vector<std::size_t>& partitionMembership();

private:
        std::size_t numElem_;
        std::vector<std::size_t> membership_;
        std::size_t maxSetIdx_;
};

#endif // __PARTITION_REFINER_HPP__
