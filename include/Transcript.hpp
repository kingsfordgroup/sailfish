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
        Sequence_(nullptr), mass_(0.0), estCount_(0.0), active_(false) { }

    Transcript(Transcript&& other) {
        id = other.id;
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        Sequence_ = other.Sequence_;
        GCCount_ = other.GCCount_;
        gcStep_ = other.gcStep_;
        gcFracLen_ = other.gcFracLen_;
        lastRegularSample_ = other.lastRegularSample_;
        mass_.store(other.mass_.load());
        estCount_.store(other.estCount_.load());
        active_ = other.active_;
    }

    Transcript& operator=(Transcript&& other) {
        id = other.id;
        RefName = std::move(other.RefName);
        RefLength = other.RefLength;
        EffectiveLength = other.EffectiveLength;
        Sequence_ = other.Sequence_;
        GCCount_ = other.GCCount_;
        gcStep_ = other.gcStep_;
        gcFracLen_ = other.gcFracLen_;
        lastRegularSample_ = other.lastRegularSample_;
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

    // Return the fractional GC content along this transcript
    // in the interval [s,e] (note; this interval is closed on both sides).
    inline int32_t gcFrac(int32_t s, int32_t e) const {
      if (gcStep_ == 1) {
        auto cs = GCCount_[s];
        auto ce = GCCount_[e];
        return std::lrint((100.0 * (ce - cs)) / (e - s + 1));
      } else {
        auto cs = gcCountInterp_(s);
        auto ce = gcCountInterp_(e);
        return std::lrint((100.0 * (ce - cs)) / (e - s + 1));
      }
    }

    void setSequence(const char* seq, bool needGC=false, uint32_t gcSampFactor=1) {
        Sequence_ = seq;
        if (needGC) { computeGCContent_(gcSampFactor); }
    }

    const char* Sequence() const { return Sequence_; }

    std::string RefName;
    uint32_t RefLength;
    double EffectiveLength;
    uint32_t id;

    double uniqueCounts{0.0};
    double totalCounts{0.0};
    double projectedCounts{0.0};
    double sharedCounts{0.0};

private:
    // NOTE: Is it worth it to check if we have GC here?
    // we should never access these without bias correction.
    inline double gcCount_(int32_t p) {
        return (gcStep_ == 1) ? static_cast<double>(GCCount_[p]) : gcCountInterp_(p);
    }
    inline double gcCount_(int32_t p) const {
        return (gcStep_ == 1) ? static_cast<double>(GCCount_[p]) : gcCountInterp_(p);
    }

    inline double gcCountInterp_(int32_t p) const {
        //std::cerr << "in gcCountInterp\n";
        if (p == RefLength - 1) {
            // If p is the last position, just return the last value
            return static_cast<double>(GCCount_.back());
        }

        // The fractional sampling factor position p would have
        double fracP = static_cast<double>(p) / gcStep_;

        // The largest sampled index for some position <= p
        uint32_t sampInd = std::floor(fracP);

        // The fraction sampling factor for the largest sampled
        // position <= p
        double fracSample = static_cast<double>(sampInd);

        int32_t nextSample{0};
        double fracNextSample{0.0};

        // special case: The last bin may not be evenly spaced.
        if (sampInd >= lastRegularSample_) {
            nextSample = GCCount_.size() - 1;
            fracNextSample = gcFracLen_;
        } else {
            nextSample = sampInd + 1;
            fracNextSample = static_cast<double>(nextSample);
        }
        double lambda = (fracP - fracSample) / (fracNextSample - fracSample);
        return lambda * GCCount_[sampInd] + (1.0 - lambda) * GCCount_[nextSample];
    }

    void computeGCContentSampled_(uint32_t step) {
        gcStep_ = step;
        const char* seq = Sequence_;
        size_t nsamp = std::ceil(static_cast<double>(RefLength) / step);
        GCCount_.reserve(nsamp + 2);

        size_t lastSamp{0};
        size_t totGC{0};
        for (size_t i = 0; i < RefLength; ++i) {
            auto c = std::toupper(seq[i]);
            if (c == 'G' or c == 'C') {
                totGC++;
            }
            if (i % step == 0) {
                GCCount_.push_back(totGC);
                lastSamp = i;
            }
        }

        if (lastSamp < RefLength - 1) {
            GCCount_.push_back(totGC);
        }

        gcFracLen_ = static_cast<double>(RefLength - 1) / gcStep_;
        lastRegularSample_ = std::ceil(gcFracLen_);
    }

    void computeGCContent_(uint32_t gcSampFactor) {
        const char* seq = Sequence_;
        GCCount_.clear();
        if (gcSampFactor == 1) {
            GCCount_.resize(RefLength, 0);
            size_t totGC{0};
            for (size_t i = 0; i < RefLength; ++i) {
                auto c = std::toupper(seq[i]);
                if (c == 'G' or c == 'C') {
                    totGC++;
                }
                GCCount_[i] = totGC;
            }
        } else {
            computeGCContentSampled_(gcSampFactor);
        }
    }

    const char* Sequence_;
    //std::unique_ptr<const char, void(*)(const char*)> Sequence =
    //    std::unique_ptr<const char, void(*)(const char*)>(nullptr, [](const char*) {});
    tbb::atomic<double> mass_;
    tbb::atomic<double> estCount_;
    bool active_;

    uint32_t gcStep_{1};
    double gcFracLen_{0.0};
    uint32_t lastRegularSample_{0};
    std::vector<uint32_t> GCCount_;
};

#endif //TRANSCRIPT
