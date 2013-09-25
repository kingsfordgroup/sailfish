#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <boost/filesystem.hpp>

namespace bfs = boost::filesystem;

int performBiasCorrection(
	bfs::path featureFile,
	bfs::path expressionFile,
	bfs::path outputFile,
	size_t numThreads);
	
void runIterativeOptimizer(int argc, char* argv[]){}

int main(int argc, char* argv[]) {
	bfs::path featureFile(argv[1]);
	bfs::path expressionFile(argv[2]);
	bfs::path outputFile(argv[3]);
	size_t numThreads = atoi(argv[4]);
	
	performBiasCorrection(featureFile, expressionFile, outputFile, numThreads);
}
