#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <unordered_map>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>

// Logger includes
#include "spdlog/spdlog.h"

#include "cuckoohash_map.hh"
#include "concurrentqueue.h"
#include "TranscriptGroup.hpp"


struct TGValue {
    TGValue(const TGValue& o) {
        weights = o.weights;
        count.store(o.count.load());
    }

    TGValue(std::vector<double>& weightIn, uint64_t countIn) :
        weights(weightIn) { count.store(countIn); }

    // const is a lie
    void normalizeAux() const {
        double sumOfAux{0.0};
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux += weights[i];
        }
        double norm = 1.0 / sumOfAux;
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] *= norm;
        }
        /* LOG SPACE
        double sumOfAux = salmon::math::LOG_0;
        for (size_t i = 0; i < weights.size(); ++i) {
            sumOfAux = salmon::math::logAdd(sumOfAux, weights[i]);
        }
        for (size_t i = 0; i < weights.size(); ++i) {
            weights[i] = std::exp(weights[i] - sumOfAux);
        }
        */
    }

    // forget synchronizing this for the time being
    mutable std::vector<double> weights;
    std::atomic<uint64_t> count{0};
};

class EquivalenceClassBuilder {
    public:
        EquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn) :
		logger_(loggerIn) {
            countMap_.reserve(1000000);
        }

        ~EquivalenceClassBuilder() {}

        void start() { active_ = true; }

        bool finish() {
            active_ = false;
            size_t totalCount{0};
            auto lt = countMap_.lock_table();
            for (auto& kv : lt) {
            //for (auto kv = countMap_.begin(); !kv.is_end(); ++kv) {
                kv.second.normalizeAux();
                totalCount += kv.second.count;
                countVec_.push_back(kv);
            }

    	    logger_->info("Computed {} rich equivalence classes "
			  "for further processing", countVec_.size());
            logger_->info("Counted {} total reads in the equivalence classes ",
                    totalCount);
            return true;
        }

        inline void insertGroup(TranscriptGroup g, uint32_t count) {
            std::vector<double> weights(g.txps.size(), 0.0);
            //auto updatefn = [count](TGValue& x) { x.count = count; };
            TGValue v(weights, count);
            countVec_.push_back(std::make_pair(g, v));
            //countMap_.upsert(g, updatefn, v);
        }

        inline void addGroup(TranscriptGroup&& g,
                             std::vector<double>& weights) {

            auto upfn = [&weights](TGValue& x) -> void {
                // update the count
                x.count++;
                // update the weights
                for (size_t i = 0; i < x.weights.size(); ++i) {
                    // Possibly atomicized in the future
                    weights[i] += x.weights[i];
                    /* LOG SPACE
                    x.weights[i] =
                        salmon::math::logAdd(x.weights[i], weights[i]);
                    */
                }
            };
            TGValue v(weights, 1);
            countMap_.upsert(g, upfn, v);
        }

        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec() {
            return countVec_;
        }

    private:
        std::atomic<bool> active_;
	    cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher> countMap_;
        std::vector<std::pair<const TranscriptGroup, TGValue>> countVec_;
    	std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_BUILDER_HPP
