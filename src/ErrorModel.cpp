#include <iostream>
#include <cstdio>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "ErrorModel.hpp"
#include "Transcript.hpp"
#include "SailfishStringUtils.hpp"
#include "UnpairedRead.hpp"
#include "ReadPair.hpp"

ErrorModel::ErrorModel(double alpha, uint32_t maxExpectedReadLen) :
    maxExpectedLen_(maxExpectedReadLen),
    mismatchLeft_(maxExpectedReadLen, AtomicMatrix<double>(16, 4, alpha)),
    mismatchRight_(maxExpectedReadLen, AtomicMatrix<double>(16, 4, alpha)),
    isEnabled_(true),
    maxLen_(0),
    burnedIn_(false) {}

bool ErrorModel::burnedIn() { return burnedIn_; }
void ErrorModel::burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

double ErrorModel::logLikelihood(bam_seq_t* read, Transcript& ref,
                                 std::vector<AtomicMatrix<double>>& mismatchProfile){
    using namespace sailfish::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    size_t transcriptIdx = bam_pos(read);
    //sailfish::stringtools::strand readStrand = (BAM_FREVERSE & read->core.flag) ? sailfish::stringtools::strand::reverse :
    //                                            sailfish::stringtools::strand::forward;
    sailfish::stringtools::strand readStrand = sailfish::stringtools::strand::forward;
    double logLike = sailfish::math::LOG_1;

    uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));//bam1_seq(read);
    auto qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));//bam1_qual(read);
    while (readIdx < bam_seq_len(read)) {
        size_t curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        size_t prevReadBase = (readIdx) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0;
        size_t refBase = samToTwoBit[ref.baseAt(transcriptIdx, readStrand)];
        size_t index = prevReadBase + refBase;
        int qval = qualStr[readIdx];
        double qual = (useQual) ? sailfish::stringtools::phredToLogProb[qval] : sailfish::math::LOG_1;
        logLike += mismatchProfile[readIdx](index, curReadBase) + qual;
        ++readIdx;
        ++transcriptIdx;
    }
    return logLike;
}

double ErrorModel::logLikelihood(const ReadPair& hit, Transcript& ref){
    double logLike = sailfish::math::LOG_1;
    if (BOOST_UNLIKELY(!isEnabled_)) { return logLike; }

    bam_seq_t* leftRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
    bam_seq_t* rightRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;

    // NOTE: Raise a warning in this case?
    if (BOOST_UNLIKELY((bam_seq_len(leftRead) > maxExpectedLen_) or
                       (bam_seq_len(rightRead) > maxExpectedLen_))) {
        return logLike;
    }

    if (leftRead) {
        logLike += logLikelihood(leftRead, ref, mismatchLeft_);
    }

    if (rightRead) {
        logLike += logLikelihood(rightRead, ref, mismatchRight_);
    }
    if (logLike == sailfish::math::LOG_0) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "error likelihood: " << logLike << "\n";
    }

    return logLike;
}

double ErrorModel::logLikelihood(const UnpairedRead& hit, Transcript& ref){
    double logLike = sailfish::math::LOG_1;
    if (BOOST_UNLIKELY(!isEnabled_)) { return logLike; }

    bam_seq_t* read = hit.read;
    // NOTE: Raise a warning in this case?
    if (BOOST_UNLIKELY(bam_seq_len(read) > maxExpectedLen_)) {
        return logLike;
    }
    logLike += logLikelihood(read, ref, mismatchLeft_);

    if (logLike == sailfish::math::LOG_0) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "error likelihood: " << logLike << "\n";
    }

    return logLike;
}

void ErrorModel::update(const UnpairedRead& hit, Transcript& ref, double p, double mass){
    if (mass == sailfish::math::LOG_0) { return; }
    if (BOOST_UNLIKELY(!isEnabled_)) { return; }
    bam_seq_t* leftRead = hit.read;
    update(leftRead, ref, p, mass, mismatchLeft_);
}

void ErrorModel::update(bam_seq_t* read, Transcript& ref, double p, double mass,
                        std::vector<AtomicMatrix<double>>& mismatchProfile) {
    using namespace sailfish::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    size_t transcriptIdx = bam_pos(read);
    //sailfish::stringtools::strand readStrand = (BAM_FREVERSE & read->core.flag) ? sailfish::stringtools::strand::reverse :
    //                                            sailfish::stringtools::strand::forward;
    sailfish::stringtools::strand readStrand = sailfish::stringtools::strand::forward;

    uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));//bam1_seq(read);
    uint8_t* qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));//bam1_qual(read);
    while (readIdx < bam_seq_len(read)) {
        size_t curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        size_t prevReadBase = (readIdx > 0) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0;
        size_t refBase = samToTwoBit[ref.baseAt(transcriptIdx, readStrand)];
        size_t index = prevReadBase + refBase;
        int qval = qualStr[readIdx];
        double qual = (useQual) ? sailfish::stringtools::phredToLogProb[qval] : sailfish::math::LOG_1;
        mismatchProfile[readIdx].increment(index, curReadBase, mass+p+qual);
        ++readIdx;
        ++transcriptIdx;
    }
    maxLen_ = std::max(maxLen_, static_cast<size_t>(bam_seq_len(read)));
    if (BOOST_UNLIKELY(maxLen_ > 250)) {
        std::lock_guard<std::mutex> lock(outputMutex_);
        std::cerr << "Encountered read longer than maximum expected length of "
                  << maxExpectedLen_ << ", not applying error model\n";
        isEnabled_ = false;
    }
}

void ErrorModel::update(const ReadPair& hit, Transcript& ref, double p, double mass){
    if (mass == sailfish::math::LOG_0) { return; }
    if (BOOST_UNLIKELY(!isEnabled_)) { return; }

    bam_seq_t* leftRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
    bam_seq_t* rightRead = (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;

    if (leftRead) {
        update(leftRead, ref, p, mass, mismatchLeft_);
    }

    if (rightRead) {
        update(rightRead, ref, p, mass, mismatchRight_);
    }
}
