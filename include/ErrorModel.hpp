#ifndef ERROR_MODEL
#define ERROR_MODEL

extern "C" {
#include "samtools/sam.h"
}

#include <mutex>
#include <atomic>
#include "tbb/concurrent_vector.h"
#include "AtomicMatrix.hpp"

struct UnpairedRead;
struct ReadPair;
class Transcript;

class ErrorModel {
public:
    ErrorModel(double alpha);
    bool burnedIn();
    void burnedIn(bool burnedIn);

    void update(const UnpairedRead&, Transcript& ref, double p, double mass);
    double logLikelihood(const UnpairedRead&, Transcript& ref);
    void update(const ReadPair&, Transcript& ref, double p, double mass);
    double logLikelihood(const ReadPair&, Transcript& ref);
private:
    void update(bam1_t* read, Transcript& ref, double p, double mass, tbb::concurrent_vector<AtomicMatrix<double>>& mismatchProfile);
    double logLikelihood(bam1_t* read, Transcript& ref, tbb::concurrent_vector<AtomicMatrix<double>>& mismatchProfile);

    tbb::concurrent_vector<AtomicMatrix<double>> mismatchLeft_;
    tbb::concurrent_vector<AtomicMatrix<double>> mismatchRight_;
    size_t maxLen_{250};
    std::atomic<bool> burnedIn_;
    std::mutex outputMutex_;
};

#endif // ERROR_MODEL
