#ifndef TRANSCRIPT
#define TRANSCRIPT

#include <atomic>
#include <cmath>
#include <limits>
#include "SailfishMath.hpp"
//#include "FragmentLengthDistribution.hpp"
#include "tbb/atomic.h"
#include "SailfishUtils.hpp"

class Transcript {
public:
    Transcript(size_t idIn, const char* name, uint32_t len) :
        RefName(name), RefLength(len), EffectiveLength(len), id(idIn),
        mass_(0.0), estCount_(0.0), active_(false) { }

    Transcript(Transcript&& other) {
        id = other.id;
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        Sequence = other.Sequence;
        mass_.store(other.mass_.load());
        estCount_.store(other.estCount_.load());
        active_ = other.active_;
    }

    Transcript& operator=(Transcript&& other) {
        id = other.id;
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        Sequence = other.Sequence;
        mass_.store(other.mass_.load());
        estCount_.store(other.estCount_.load());
        active_ = other.active_;
        return *this;
    }

    inline double estCount() { return estCount_.load(); }
    //inline size_t uniqueCount() { return uniqueCount_.load(); }
    //inline size_t totalCount() { return totalCount_.load(); }
    //inline void addUniqueCount(size_t newCount) { uniqueCount_ += newCount; }
    //inline void addTotalCount(size_t newCount) { totalCount_ += newCount; }

    inline void setEstCount(double sc) {
        estCount_.store(sc);
    }

    inline void addEstCount(double sc) {
	    sailfish::utils::incLoop(estCount_, sc);
    }

    inline void addMass(double mass) {
	    sailfish::utils::incLoop(mass_, mass);
    }

    inline void setMass(double mass) {
        mass_.store(mass);
    }

    inline double mass(bool fauxArg) {
        return mass_.load();
    }


    inline double mass() {
        return mass_.load();
    }

    void setActive() { active_ = true; }
    bool getActive() { return active_; }
    /**
      *  NOTE: Adopted from "est_effective_length" at (https://github.com/adarob/eXpress/blob/master/src/targets.cpp)
      *  originally written by Adam Roberts.
      *
      *
      */
    /*
    double updateEffectiveLength(const FragmentLengthDistribution& fragLengthDist) {

        double effectiveLength = salmon::math::LOG_0;
        double refLen = static_cast<double>(RefLength);
        double logLength = std::log(refLen);

        if (logLength < fragLengthDist.mean()) {
            effectiveLength = logLength;
        } else {
            uint32_t mval = fragLengthDist.maxVal();
            for (size_t l = fragLengthDist.minVal(); l <= std::min(RefLength, mval); ++l) {
                effectiveLength = salmon::math::logAdd(
                        effectiveLength,
                        fragLengthDist.pmf(l) + std::log(refLen - l + 1));
            }
        }

        return effectiveLength;
    }

    double getCachedEffectiveLength() {
        return cachedEffectiveLength_.load();
    }

    double getEffectiveLength(const FragmentLengthDistribution& fragLengthDist,
                              size_t currObs,
                              size_t burnInObs) {
        if (lastUpdate_ == 0 or
            (currObs - lastUpdate_ >= 250000) or
            (lastUpdate_ < burnInObs and currObs > burnInObs)) {
            // compute new number
            double cel = updateEffectiveLength(fragLengthDist);
            cachedEffectiveLength_.store(cel);
            lastUpdate_.store(currObs);
            //priorMass_ = cel + logPerBasePrior_;
            return cachedEffectiveLength_.load();
        } else {
            // return cached number
            return cachedEffectiveLength_.load();
        }
    }
    */

    std::string RefName;
    uint32_t RefLength;
    double EffectiveLength;
    uint32_t id;

    double uniqueCounts{0.0};
    double totalCounts{0.0};
    double projectedCounts{0.0};
    double sharedCounts{0.0};
    std::string Sequence;

private:
    tbb::atomic<double> mass_;
    tbb::atomic<double> estCount_;
    bool active_;
};

#endif //TRANSCRIPT
