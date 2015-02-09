#ifndef __FORGETTING_MASS_CALCULATOR__
#define __FORGETTING_MASS_CALCULATOR__

#include "SpinLock.hpp"
#include "SailfishMath.hpp"

class ForgettingMassCalculator {
    public:
    ForgettingMassCalculator(double forgettingFactor = 0.65) :
        batchNum_(0), forgettingFactor_(forgettingFactor),
        logForgettingMass_(sailfish::math::LOG_1) {}

    ForgettingMassCalculator(const ForgettingMassCalculator&) = delete;
    ForgettingMassCalculator(ForgettingMassCalculator&&) = default;
    ForgettingMassCalculator& operator=(ForgettingMassCalculator&&) = default;
    ForgettingMassCalculator& operator=(const ForgettingMassCalculator&) = delete;

    double operator()() {
#if defined __APPLE__
        spin_lock::scoped_lock sl(ffMutex_);
#else
        std::lock_guard<std::mutex> lock(ffMutex_);
#endif
        ++batchNum_;
        if (batchNum_ > 1) {
            logForgettingMass_ += forgettingFactor_ * std::log(static_cast<double>(batchNum_-1)) -
                std::log(std::pow(static_cast<double>(batchNum_), forgettingFactor_) - 1);
        }
        return logForgettingMass_;
    }

    private:
        uint64_t batchNum_;
        double forgettingFactor_;
        double logForgettingMass_;
#if defined __APPLE__
        spin_lock ffMutex_;
#else
        std::mutex ffMutex_;
#endif
};

#endif //__FORGETTING_MASS_CALCULATOR__
