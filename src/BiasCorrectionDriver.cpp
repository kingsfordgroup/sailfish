#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

int performBiasCorrection(
        bfs::path featureFile,
        bfs::path expressionFile,
        double estimatedReadLength,
        double kmersPerRead,
        uint64_t mappedKmers,
        uint32_t merLen,
        bfs::path outputFile,
        size_t numThreads);

void runIterativeOptimizer(int argc, char* argv[]){}

int main(int argc, char* argv[]) {
        bfs::path featureFile(argv[1]);
        bfs::path expressionFile(argv[2]);
        double estimatedReadLength = atod(argv[3]);
        double kmersPerRead = atod(argv[4]);
        uint64_t mappedKmers = atol(argv[5]);
        uint32_t mappedKmers = atoi(argv[6]);
        bfs::path outputFile(argv[7]);
        size_t numThreads = atoi(argv[8]);

        performBiasCorrection(featureFile, expressionFile, estimatedReadLength, kmersPerRead,
                              mappedKmers, merLen, outputFile, numThreads);
}
