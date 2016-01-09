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

    uint32_t gcCount(int32_t p) { return GCCount_[p]; }
    uint32_t gcCount(int32_t p) const { return GCCount_[p]; }

    void setSequence(const char* seq) {
        Sequence_ = seq;
        GCCount_.clear();
        GCCount_.resize(RefLength, 0);
        size_t totGC{0};
        for (size_t i = 0; i < RefLength; ++i) {
            auto c = std::toupper(seq[i]);
            if (c == 'G' or c == 'C') {
                totGC++;
            }
            GCCount_[i] = totGC;
        }
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
    const char* Sequence_;
    std::vector<uint32_t> GCCount_;
    //std::unique_ptr<const char, void(*)(const char*)> Sequence =
    //    std::unique_ptr<const char, void(*)(const char*)>(nullptr, [](const char*) {});
    tbb::atomic<double> mass_;
    tbb::atomic<double> estCount_;
    bool active_;
};

#endif //TRANSCRIPT
