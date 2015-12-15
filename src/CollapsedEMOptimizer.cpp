#include <vector>
#include <unordered_map>
#include <atomic>

#include "tbb/task_scheduler_init.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/partitioner.h"
#include "concurrentqueue.h"

#include <boost/math/special_functions/digamma.hpp>

// C++ string formatting library
#include "spdlog/details/format.h"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

#include "CollapsedEMOptimizer.hpp"
#include "Transcript.hpp"
#include "TranscriptGroup.hpp"
#include "SailfishMath.hpp"
#include "ReadExperiment.hpp"
#include "BootstrapWriter.hpp"
#include "MultinomialSampler.hpp"

using BlockedIndexRange =  tbb::blocked_range<size_t>;

// intelligently chosen value adopted from
// https://github.com/pachterlab/kallisto/blob/master/src/EMAlgorithm.h#L18
constexpr double minEQClassWeight = std::numeric_limits<double>::denorm_min();
constexpr double minWeight = std::numeric_limits<double>::denorm_min();

template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff) {
    double alphaSum = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i) {
        if (alphas[i] <= cutoff) { alphas[i] = 0.0; }
        alphaSum += alphas[i];
    }
    return alphaSum;
}

double normalize(std::vector<tbb::atomic<double>>& vec) {
    double sum{0.0};
    for (auto& v : vec) {
        sum += v;
    }

    // too small!
    if (sum < minWeight) {
        return sum;
    }

    double invSum  = 1.0 / sum;
    for (auto& v : vec) {
        v.store(v.load() * invSum);
    }

    return sum;
}

/*
 * Use atomic compare-and-swap to update val to
 * val + inc.  Update occurs in a loop in case other
 * threads update in the meantime.
 */
void incLoop(tbb::atomic<double>& val, double inc) {
        double oldMass = val.load();
        double returnedMass = oldMass;
        double newMass{oldMass + inc};
        do {
            oldMass = returnedMass;
            newMass = oldMass + inc;
            returnedMass = val.compare_and_swap(newMass, oldMass);
        } while (returnedMass != oldMass);
}

/*
 * Same as above, but overloaded for "plain" doubles
 *
 */
void incLoop(double& val, double inc) {
    val += inc;
}

/**
 * Single-threaded EM-update routine for use in bootstrapping
 */
template <typename VecT>
void EMUpdate_(
        std::vector<std::vector<uint32_t>>& txpGroupLabels,
        std::vector<std::vector<double>>& txpGroupWeights,
        std::vector<uint64_t>& txpGroupCounts,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const VecT& alphaIn,
        VecT& alphaOut) {

    assert(alphaIn.size() == alphaOut.size());

    size_t numEqClasses = txpGroupLabels.size();
    for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
        uint64_t count = txpGroupCounts[eqID];
        // for each transcript in this class
        const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
        const std::vector<double>& auxs = txpGroupWeights[eqID];

        double denom = 0.0;
        size_t groupSize = txps.size();
        // If this is a single-transcript group,
        // then it gets the full count.  Otherwise,
        // update according to our VBEM rule.
        if (BOOST_LIKELY(groupSize > 1)) {
           for (size_t i = 0; i < groupSize; ++i) {
               auto tid = txps[i];
               auto aux = auxs[i];
               double v = alphaIn[tid] * aux;
               denom += v;
            }

            if (denom <= ::minEQClassWeight) {
                // tgroup.setValid(false);
            } else {
                double invDenom = count / denom;
                for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    double v = alphaIn[tid] * aux;
                    if (!std::isnan(v)) {
                        incLoop(alphaOut[tid], v * invDenom);
                    }
                }
            }
        } else {
            incLoop(alphaOut[txps.front()], count);
        }
    }
}

/**
 * Single-threaded VBEM-update routine for use in bootstrapping
 */
template <typename VecT>
void VBEMUpdate_(
	std::vector<std::vector<uint32_t>>& txpGroupLabels,
	std::vector<std::vector<double>>& txpGroupWeights,
	std::vector<uint64_t>& txpGroupCounts,
	std::vector<Transcript>& transcripts,
	Eigen::VectorXd& effLens,
        double priorAlpha,
        double totLen,
        const VecT& alphaIn,
        VecT& alphaOut,
	VecT& expTheta) {

    assert(alphaIn.size() == alphaOut.size());

    size_t numEQClasses = txpGroupLabels.size();
    double alphaSum = {0.0};
    for (auto& e : alphaIn) { alphaSum += e; }

    double logNorm = boost::math::digamma(alphaSum);


    double prior = priorAlpha;
    double priorNorm = prior * totLen;

    for (size_t i = 0; i < transcripts.size(); ++i) {
	if (alphaIn[i] > ::minWeight) {
	    expTheta[i] = std::exp(boost::math::digamma(alphaIn[i]) - logNorm);
	} else {
	    expTheta[i] = 0.0;
	}
	alphaOut[i] = prior;
    }

    for (size_t eqID = 0; eqID < numEQClasses; ++eqID) {
	uint64_t count = txpGroupCounts[eqID];
	const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
	const std::vector<double>& auxs = txpGroupWeights[eqID];

	double denom = 0.0;
	size_t groupSize = txps.size();
	// If this is a single-transcript group,
	// then it gets the full count.  Otherwise,
	// update according to our VBEM rule.
	if (BOOST_LIKELY(groupSize > 1)) {
	    for (size_t i = 0; i < groupSize; ++i) {
		auto tid = txps[i];
		auto aux = auxs[i];
		if (expTheta[tid] > 0.0) {
		    double v = expTheta[tid] * aux;
		    denom += v;
		}
	    }
	    if (denom <= ::minEQClassWeight) {
		// tgroup.setValid(false);
	    } else {
		double invDenom = count / denom;
		for (size_t i = 0; i < groupSize; ++i) {
		    auto tid = txps[i];
		    auto aux = auxs[i];
		    if (expTheta[tid] > 0.0) {
			double v = expTheta[tid] * aux;
			incLoop(alphaOut[tid], v * invDenom);
		    }
		}
	    }

	} else {
	    incLoop(alphaOut[txps.front()], count);
	}
    }
}

/*
 * Use the "standard" EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
void EMUpdate_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const CollapsedEMOptimizer::VecType& alphaIn,
        CollapsedEMOptimizer::VecType& alphaOut) {

    assert(alphaIn.size() == alphaOut.size());

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &alphaIn, &transcripts, &alphaOut](const BlockedIndexRange& range) -> void {
            for (auto eqID : boost::irange(range.begin(), range.end())) {
            auto& kv = eqVec[eqID];

            uint64_t count = kv.second.count;
            // for each transcript in this class
            const TranscriptGroup& tgroup = kv.first;
            if (tgroup.valid) {
                const std::vector<uint32_t>& txps = tgroup.txps;
                const std::vector<double>& auxs = kv.second.weights;

                double denom = 0.0;
                size_t groupSize = txps.size();
                // If this is a single-transcript group,
                // then it gets the full count.  Otherwise,
                // update according to our VBEM rule.
                if (BOOST_LIKELY(groupSize > 1)) {
                    for (size_t i = 0; i < groupSize; ++i) {
                    auto tid = txps[i];
                    auto aux = auxs[i];
                    //double el = effLens(tid);
                    //if (el <= 0) { el = 1.0; }
                    double v = alphaIn[tid] * aux;
                    denom += v;
                    }

                    if (denom <= ::minEQClassWeight) {
                        // tgroup.setValid(false);
                    } else {

                        double invDenom = count / denom;
                        for (size_t i = 0; i < groupSize; ++i) {
                            auto tid = txps[i];
                            auto aux = auxs[i];
                            double v = alphaIn[tid] * aux;
                            if (!std::isnan(v)) {
                                incLoop(alphaOut[tid], v * invDenom);
                            }
                        }
                    }
                } else {
                    incLoop(alphaOut[txps.front()], count);
                }
            }
    }
    });

}

/*
 * Use the Variational Bayesian EM algorithm over equivalence
 * classes to estimate the latent variables (alphaOut)
 * given the current estimates (alphaIn).
 */
void VBEMUpdate_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        double priorAlpha,
        double totLen,
        const CollapsedEMOptimizer::VecType& alphaIn,
        CollapsedEMOptimizer::VecType& alphaOut,
	    CollapsedEMOptimizer::VecType& expTheta) {

    assert(alphaIn.size() == alphaOut.size());

    double alphaSum = {0.0};
    for (auto& e : alphaIn) { alphaSum += e; }

    double logNorm = boost::math::digamma(alphaSum);

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(transcripts.size())),
            [logNorm, priorAlpha, totLen, &effLens, &alphaIn,
             &alphaOut, &expTheta]( const BlockedIndexRange& range) -> void {

             double prior = priorAlpha;
             double priorNorm = prior * totLen;

             for (auto i : boost::irange(range.begin(), range.end())) {
                if (alphaIn[i] > ::minWeight) {
                    expTheta[i] = std::exp(boost::math::digamma(alphaIn[i].load()) - logNorm);
                } else {
                    expTheta[i] = 0.0;
                }
                alphaOut[i] = prior;
            }
        });

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &alphaIn,
             &alphaOut, &expTheta]( const BlockedIndexRange& range) -> void {
            for (auto eqID : boost::irange(range.begin(), range.end())) {
            auto& kv = eqVec[eqID];

            uint64_t count = kv.second.count;
            // for each transcript in this class
            const TranscriptGroup& tgroup = kv.first;
            if (tgroup.valid) {
                const std::vector<uint32_t>& txps = tgroup.txps;
                const std::vector<double>& auxs = kv.second.weights;

                double denom = 0.0;
                size_t groupSize = txps.size();
                // If this is a single-transcript group,
                // then it gets the full count.  Otherwise,
                // update according to our VBEM rule.
                if (BOOST_LIKELY(groupSize > 1)) {
                    for (size_t i = 0; i < groupSize; ++i) {
                        auto tid = txps[i];
                        auto aux = auxs[i];
                        if (expTheta[tid] > 0.0) {
                            double v = expTheta[tid] * aux;
                            denom += v;
                       }
                    }
                    if (denom <= ::minEQClassWeight) {
                        // tgroup.setValid(false);
                    } else {
                        double invDenom = count / denom;
                        for (size_t i = 0; i < groupSize; ++i) {
                            auto tid = txps[i];
                            auto aux = auxs[i];
                            if (expTheta[tid] > 0.0) {
                              double v = expTheta[tid] * aux;
                              incLoop(alphaOut[tid], v * invDenom);
                            }
                        }
                    }

                } else {
                    incLoop(alphaOut[txps.front()], count);
                }
            }
        }});

}

template <typename VecT>
size_t markDegenerateClasses(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        VecT& alphaIn,
        std::shared_ptr<spdlog::logger> jointLog,
        bool verbose=false) {

    size_t numDropped{0};
    size_t idx{0};
    for (auto& kv : eqVec) {
        uint64_t count = kv.second.count;
        // for each transcript in this class
        const TranscriptGroup& tgroup = kv.first;
        const std::vector<uint32_t>& txps = tgroup.txps;
        const std::vector<double>& auxs = kv.second.weights;

        double denom = 0.0;
        for (size_t i = 0; i < txps.size(); ++i) {
            auto tid = txps[i];
            auto aux = auxs[i];
            double v = alphaIn[tid] * aux;
            if (!std::isnan(v)) {
                denom += v;
            } else {
                std::cerr << "val is NAN; alpha( "
                          << tid << " ) = " << alphaIn[tid]
                          << ", aux = " << aux << "\n";
            }
        }
        if (denom <= minEQClassWeight) {
            fmt::MemoryWriter errstream;

            errstream << "\nDropping weighted eq class\n";
            errstream << "============================\n";

            errstream << "denom = 0, count = " << count << "\n";
            errstream << "class = { ";
            for (auto e : txps) {
                errstream << e << " ";
            }
            errstream << "}\n";
            errstream << "alphas = { ";
            for (auto e : txps) {
                errstream << alphaIn[e] << " ";
            }
            errstream << "}\n";
            errstream << "weights = { ";
            for (auto e : auxs) {
                errstream << e << " ";
            }
            errstream << "}\n";
            errstream << "============================\n\n";

            bool verbose{false};
            if (verbose) {
                jointLog->info(errstream.str());
            }
            ++numDropped;
            kv.first.setValid(false);
        }
    }
    return numDropped;
}


CollapsedEMOptimizer::CollapsedEMOptimizer() {}

bool doBootstrap(
        std::vector<std::vector<uint32_t>>& txpGroups,
        std::vector<std::vector<double>>& txpGroupWeights,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        std::vector<double>& sampleWeights,
        uint64_t totalNumFrags,
        double uniformTxpWeight,
        std::atomic<uint32_t>& bsNum,
        SailfishOpts& sopt,
        BootstrapWriter* bootstrapWriter,
        double relDiffTolerance,
        uint32_t maxIter) {


    auto& jointLog = sopt.jointLog;
    bool useVBEM{sopt.useVBOpt};
    size_t numClasses = txpGroups.size();
    CollapsedEMOptimizer::SerialVecType alphas(transcripts.size(), 0.0);
    CollapsedEMOptimizer::SerialVecType alphasPrime(transcripts.size(), 0.0);
    CollapsedEMOptimizer::SerialVecType expTheta(transcripts.size(), 0.0);
    std::vector<uint64_t> sampCounts(numClasses, 0);

    uint32_t numBootstraps = sopt.numBootstraps;

    std::random_device rd;
    MultinomialSampler msamp(rd);

    while (bsNum++ < numBootstraps) {
        // Do a new bootstrap
        msamp(sampCounts.begin(), totalNumFrags, numClasses, sampleWeights.begin());

        double totalLen{0.0};
        for (size_t i = 0; i < transcripts.size(); ++i) {
            alphas[i] = transcripts[i].getActive() ? uniformTxpWeight * totalNumFrags : 0.0;
            totalLen += effLens(i);
        }

        bool converged{false};
        double maxRelDiff = -std::numeric_limits<double>::max();
        size_t itNum = 0;

        // If we use VBEM, we'll need the prior parameters
        double priorAlpha = 0.01;
        double minAlpha = 1e-8;
        double alphaCheckCutoff = 1e-2;
        double cutoff = (useVBEM) ? (priorAlpha + minAlpha) : minAlpha;

        while (itNum < maxIter and !converged) {

            if (useVBEM) {
                VBEMUpdate_(txpGroups, txpGroupWeights, sampCounts, transcripts,
                        effLens, priorAlpha, totalLen, alphas, alphasPrime, expTheta);
            } else {
                EMUpdate_(txpGroups, txpGroupWeights, sampCounts, transcripts,
                        effLens, alphas, alphasPrime);
            }

            converged = true;
            maxRelDiff = -std::numeric_limits<double>::max();
            for (size_t i = 0; i < transcripts.size(); ++i) {
                if (alphas[i] > alphaCheckCutoff) {
                    double relDiff = std::abs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
                    maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
                    if (relDiff > relDiffTolerance) {
                        converged = false;
                    }
                }
                alphas[i] = alphasPrime[i];
                alphasPrime[i] = 0.0;
            }

            ++itNum;
        }

        // Truncate tiny expression values
        double alphaSum = truncateCountVector(alphas, cutoff);

        if (alphaSum < minWeight) {
            jointLog->error("Total alpha weight was too small! "
                    "Make sure you ran sailfish correctly.");
            return false;
        }

        bootstrapWriter->writeBootstrap(alphas);
    }
    return true;
}

void updateEqClassWeights(std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
                          Eigen::VectorXd& effLens) {
    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &effLens]( const BlockedIndexRange& range) -> void {
                // For each index in the equivalence class vector
                for (auto eqID : boost::irange(range.begin(), range.end())) {
                    // The vector entry
                    auto& kv = eqVec[eqID];
                    // The label of the equivalence class
                    const TranscriptGroup& k = kv.first;
                    // The size of the label
                    size_t classSize = k.txps.size();
                    // The weights of the label
                    TGValue& v = kv.second;

                    // Iterate over each weight and set it equal to
                    // 1 / effLen of the corresponding transcript
                    double wsum{0.0};
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] = (kv.second.count / effLens(k.txps[i]));
                        wsum += v.weights[i];
                    }
                    double wnorm = 1.0 / wsum;
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] *= wnorm;
                    }
                }
            });
}

bool CollapsedEMOptimizer::gatherBootstraps(
        ReadExperiment& readExp,
        SailfishOpts& sopt,
        BootstrapWriter* bootstrapWriter,
        double relDiffTolerance,
        uint32_t maxIter) {

    std::vector<Transcript>& transcripts = readExp.transcripts();
    using VecT = CollapsedEMOptimizer::SerialVecType;

    VecT alphas(transcripts.size(), 0.0);
    VecT alphasPrime(transcripts.size(), 0.0);
    VecT expTheta(transcripts.size());
    Eigen::VectorXd effLens(transcripts.size());

    uint32_t numBootstraps = sopt.numBootstraps;

    // Fill in the effective length vector
    double totalLen{0.0};
    for (size_t i = 0; i < transcripts.size(); ++i) {
        effLens(i) = (sopt.noEffectiveLengthCorrection) ?
                        transcripts[i].RefLength : transcripts[i].EffectiveLength;
        if (effLens(i) <= 1.0) { effLens(i) = 1.0; }
        totalLen += effLens(i);
    }

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    std::unordered_set<uint32_t> activeTranscriptIDs;
    for (auto& kv : eqVec) {
        auto& tg = kv.first;
        for (auto& t : tg.txps) {
            transcripts[t].setActive();
            activeTranscriptIDs.insert(t);
        }
    }

    bool useVBEM{sopt.useVBOpt};
    // If we use VBEM, we'll need the prior parameters
    double priorAlpha = 0.01;

    auto jointLog = sopt.jointLog;

    jointLog->info("Will draw {} bootstrap samples", numBootstraps);
    jointLog->info("Optimizing over {} equivalence classes", eqVec.size());

    double totalNumFrags{static_cast<double>(readExp.numMappedFragments())};

    if (activeTranscriptIDs.size() == 0) {
        jointLog->error("It seems that no transcripts are expressed; something is likely wrong!");
        jointLog->flush();
        return false;
    }

    double scale = 1.0 / activeTranscriptIDs.size();
    for (size_t i = 0; i < transcripts.size(); ++i) {
        alphas[i] = transcripts[i].getActive() ? scale * totalNumFrags : 0.0;
    }

    auto numRemoved = markDegenerateClasses(eqVec, alphas, sopt.jointLog);
    sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
            numRemoved);

    size_t itNum{0};
    double minAlpha = 1e-8;
    double cutoff = (useVBEM) ? (priorAlpha + minAlpha) : minAlpha;

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &effLens]( const BlockedIndexRange& range) -> void {
                // For each index in the equivalence class vector
                for (auto eqID : boost::irange(range.begin(), range.end())) {
                    // The vector entry
                    auto& kv = eqVec[eqID];
                    // The label of the equivalence class
                    const TranscriptGroup& k = kv.first;
                    // The size of the label
                    size_t classSize = k.txps.size();
                    // The weights of the label
                    TGValue& v = kv.second;

                    // Iterate over each weight and set it equal to
                    // 1 / effLen of the corresponding transcript
                    double wsum{0.0};
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] = (kv.second.count / effLens(k.txps[i]));
                        wsum += v.weights[i];
                    }
                    double wnorm = 1.0 / wsum;
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] *= wnorm;
                    }
                }
            });

    // Since we will use the same weights and transcript groups for each
    // of the bootstrap samples (only the count vector will change), it
    // makes sense to keep only one copy of these.
    using TGroupLabelT = std::vector<uint32_t>;
    using TGroupWeightVec = std::vector<double>;
    std::vector<TGroupLabelT> txpGroups;
    std::vector<TGroupWeightVec> txpGroupWeights;
    std::vector<uint64_t> origCounts;
    uint64_t totalCount{0};

    for (auto& kv : eqVec) {
        uint64_t count = kv.second.count;
        // for each transcript in this class
        const TranscriptGroup& tgroup = kv.first;
        if (tgroup.valid) {
            const std::vector<uint32_t>& txps = tgroup.txps;
            const std::vector<double>& auxs = kv.second.weights;
            txpGroups.push_back(txps);
            txpGroupWeights.push_back(auxs);
            origCounts.push_back(count);
            totalCount += count;
        }
    }

    double floatCount = totalCount;
    std::vector<double> samplingWeights(txpGroups.size(), 0.0);
    for (size_t i = 0; i < origCounts.size(); ++i) {
        samplingWeights[i] = origCounts[i] / floatCount;
    }

    size_t numWorkerThreads{1};
    if (sopt.numThreads > 1 and numBootstraps > 1) {
        numWorkerThreads = std::min(sopt.numThreads - 1, numBootstraps - 1);
    }

    std::atomic<uint32_t> bsCounter{0};
    std::vector<std::thread> workerThreads;
    for (size_t tn = 0; tn < numWorkerThreads; ++tn) {
        workerThreads.emplace_back(doBootstrap,
                std::ref(txpGroups),
                std::ref(txpGroupWeights),
                std::ref(transcripts),
                std::ref(effLens),
                std::ref(samplingWeights),
                totalCount,
                scale,
                std::ref(bsCounter),
                std::ref(sopt),
                bootstrapWriter,
                relDiffTolerance,
                maxIter);
    }

    for (auto& t : workerThreads) {
        t.join();
    }
    return true;
}

bool CollapsedEMOptimizer::optimize(ReadExperiment& readExp,
        SailfishOpts& sopt,
        double relDiffTolerance,
        uint32_t maxIter) {

    uint32_t minIter = 50;
    bool doBiasCorrect = sopt.biasCorrect;

    tbb::task_scheduler_init tbbScheduler(sopt.numThreads);
    std::vector<Transcript>& transcripts = readExp.transcripts();

    using VecT = CollapsedEMOptimizer::VecType;
    // With atomics
    VecType alphas(transcripts.size(), 0.0);
    VecType alphasPrime(transcripts.size(), 0.0);
    VecType expTheta(transcripts.size());
    VecType uniqueCounts(transcripts.size(), 0.0);
    Eigen::VectorXd effLens(transcripts.size());

    // Fill in the effective length vector
    double totalLen{0.0};
    for (size_t i = 0; i < transcripts.size(); ++i) {
        effLens(i) = (sopt.noEffectiveLengthCorrection) ?
                        transcripts[i].RefLength : transcripts[i].EffectiveLength;
        if (effLens(i) <= 1.0) { effLens(i) = 1.0; }
        totalLen += effLens(i);
    }

    std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec =
        readExp.equivalenceClassBuilder().eqVec();

    tbb::parallel_for(BlockedIndexRange(size_t(0), size_t(eqVec.size())),
            [&eqVec, &effLens, &transcripts]( const BlockedIndexRange& range) -> void {
                // For each index in the equivalence class vector
                for (auto eqID : boost::irange(range.begin(), range.end())) {
                    // The vector entry
                    auto& kv = eqVec[eqID];
                    // The label of the equivalence class
                    const TranscriptGroup& k = kv.first;
                    // The size of the label
                    size_t classSize = k.txps.size();
                    // The weights of the label
                    TGValue& v = kv.second;

                    // Iterate over each weight and set it equal to
                    // 1 / effLen of the corresponding transcript
                    double wsum{0.0};
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] = (kv.second.count / effLens(k.txps[i]));
                        wsum += v.weights[i];
                    }

                    double wnorm = 1.0 / wsum;
                    for (size_t i = 0; i < classSize; ++i) {
                        v.weights[i] *= wnorm;
                    }

                }
            });

    std::unordered_set<uint32_t> activeTranscriptIDs;
    for (auto& kv : eqVec) {
        auto& tg = kv.first;
        for (size_t i = 0; i < tg.txps.size(); ++i) {
            auto& t = tg.txps[i];
            transcripts[t].setActive();
            activeTranscriptIDs.insert(t);
        }
    }

    bool useVBEM{sopt.useVBOpt};
    // If we use VBEM, we'll need the prior parameters
    double priorAlpha = 0.01;

    auto jointLog = sopt.jointLog;

    jointLog->info("Optimizing over {} equivalence classes", eqVec.size());

    double totalNumFrags{static_cast<double>(readExp.numMappedFragments())};

    if (activeTranscriptIDs.size() == 0) {
        jointLog->error("It seems that no transcripts are expressed; something is likely wrong!");
        jointLog->flush();
        return false;
    }

    double scale = 1.0 / activeTranscriptIDs.size();
    for (size_t i = 0; i < transcripts.size(); ++i) {
        alphas[i] = transcripts[i].getActive() ? scale * totalNumFrags : 0.0;
    }

    //auto numRemoved = markDegenerateClasses(eqVec, alphas, sopt.jointLog);
    //sopt.jointLog->info("Marked {} weighted equivalence classes as degenerate",
    //        numRemoved);

    size_t itNum{0};
    double minAlpha = 1e-8;
    double alphaCheckCutoff = 1e-2;
    double cutoff = (useVBEM) ? (priorAlpha + minAlpha) : minAlpha;

    bool converged{false};
    double maxRelDiff = -std::numeric_limits<double>::max();
    while (itNum < minIter or (itNum < maxIter and !converged)) {

        // Recompute the effective lengths to account for sequence-specific
        // bias.  Consider a better metric here.
        if (doBiasCorrect and (itNum == minIter or ((itNum + minIter) % 500) == 0)) {
            jointLog->info("iteration {}, recomputing effective lengths", itNum);
            effLens = sailfish::utils::updateEffectiveLengths(
                        readExp,
                        effLens,
                        alphas);
            // Check for strangeness with the lengths.
            for (size_t i = 0; i < effLens.size(); ++i) {
                if (effLens(i) <= 0.0) {
                    jointLog->warn("Transcript {} had length {}", i, effLens(i));
                }
            }
            updateEqClassWeights(eqVec, effLens);
        }

        if (useVBEM) {
            VBEMUpdate_(eqVec, transcripts, effLens,
                        priorAlpha, totalLen, alphas, alphasPrime, expTheta);
        } else {
            EMUpdate_(eqVec, transcripts, effLens, alphas, alphasPrime);
        }

        converged = true;
        maxRelDiff = -std::numeric_limits<double>::max();
        for (size_t i = 0; i < transcripts.size(); ++i) {
            if (alphasPrime[i] > alphaCheckCutoff) {
                double relDiff = std::fabs(alphas[i] - alphasPrime[i]) / alphasPrime[i];
                maxRelDiff = (relDiff > maxRelDiff) ? relDiff : maxRelDiff;
                if (relDiff > relDiffTolerance) {
                    converged = false;
                }
            }
            alphas[i] = alphasPrime[i];
            alphasPrime[i] = 0.0;
        }

        if (itNum % 100 == 0) {
            jointLog->info("iteration = {} | max rel diff. = {}",
                            itNum, maxRelDiff);
        }

        ++itNum;
    }

    jointLog->info("iteration = {} | max rel diff. = {}",
                    itNum, maxRelDiff);

    // Truncate tiny expression values
    double alphaSum = truncateCountVector(alphas, cutoff);

    if (alphaSum < minWeight) {
        jointLog->error("Total alpha weight was too small! "
                        "Make sure you ran sailfish correctly.");
        return false;
    }

    // Set the mass of each transcript using the
    // computed alphas.
    for (size_t i = 0; i < transcripts.size(); ++i) {
        // Set the mass to the normalized (after truncation)
        // relative abundance
        if (doBiasCorrect) { transcripts[i].EffectiveLength = effLens(i); }
        transcripts[i].setEstCount(alphas[i]);
        transcripts[i].setMass(alphas[i] / alphaSum);
    }
    return true;
}
