#include <iostream>
#include <cstdio>
#include "ErrorModel.hpp"
#include "Transcript.hpp"
#include "SailfishStringUtils.hpp"
#include "UnpairedRead.hpp"
#include "ReadPair.hpp"

ErrorModel::ErrorModel(double alpha) :
    mismatchLeft_(250, AtomicMatrix<double>(16, 4, alpha)),
    mismatchRight_(250, AtomicMatrix<double>(16, 4, alpha)),
    maxLen_(0),
    burnedIn_(false) {}

bool ErrorModel::burnedIn() { return burnedIn_; }
void ErrorModel::burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

double ErrorModel::logLikelihood(bam1_t* read, Transcript& ref,
                                 tbb::concurrent_vector<AtomicMatrix<double>>& mismatchProfile){
    using namespace sailfish::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    size_t transcriptIdx = read->core.pos;
    //sailfish::stringtools::strand readStrand = (bam1_strand(read)) ? sailfish::stringtools::strand::reverse :
    //    sailfish::stringtools::strand::forward;
    sailfish::stringtools::strand readStrand = sailfish::stringtools::strand::forward;
    double logLike = sailfish::math::LOG_1;

    uint8_t* qseq = bam_get_seq(read);//bam1_seq(read);
    auto qualStr = bam_get_qual(read);//bam1_qual(read);
    while (readIdx < read->core.l_qseq) {
        size_t curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];// samToTwoBit[bam1_seqi(qseq, readIdx)];
        size_t prevReadBase = (readIdx) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0; //(readIdx) ? samToTwoBit[bam1_seqi(qseq, readIdx-1)] : 0;
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
    bam1_t* leftRead = (hit.read1->core.mpos < hit.read2->core.mpos) ? hit.read1 : hit.read2;
    bam1_t* rightRead = (hit.read1->core.mpos < hit.read2->core.mpos) ? hit.read2 : hit.read1;

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
    bam1_t* read = hit.read;
    logLike += logLikelihood(read, ref, mismatchLeft_);

    if (logLike == sailfish::math::LOG_0) {
            std::lock_guard<std::mutex> lock(outputMutex_);
            std::cerr << "error likelihood: " << logLike << "\n";
    }

    return logLike;
}

void ErrorModel::update(const UnpairedRead& hit, Transcript& ref, double p, double mass){
    if (mass == sailfish::math::LOG_0) { return; }
    bam1_t* leftRead = hit.read;
    update(leftRead, ref, p, mass, mismatchLeft_);
}

void ErrorModel::update(bam1_t* read, Transcript& ref, double p, double mass,
                        tbb::concurrent_vector<AtomicMatrix<double>>& mismatchProfile) {
    using namespace sailfish::stringtools;
    bool useQual{false};
    size_t readIdx{0};
    size_t transcriptIdx = read->core.pos;
    //sailfish::stringtools::strand readStrand = (bam1_strand(read)) ? sailfish::stringtools::strand::reverse :
    //    sailfish::stringtools::strand::forward;
    sailfish::stringtools::strand readStrand = sailfish::stringtools::strand::forward;

    uint8_t* qseq = bam_get_seq(read);//bam1_seq(read);
    uint8_t* qualStr = bam_get_qual(read);//bam1_qual(read);
    while (readIdx < read->core.l_qseq) {
        size_t curReadBase = bam_seqi(qseq, readIdx);//samToTwoBit[bam1_seqi(qseq, readIdx)];
        size_t prevReadBase = (readIdx > 0) ? samToTwoBit[bam_seqi(qseq, readIdx-1)] : 0;//(readIdx > 0) ? samToTwoBit[bam1_seqi(qseq, readIdx-1)] : 0;
        size_t refBase = samToTwoBit[ref.baseAt(transcriptIdx, readStrand)];
        size_t index = prevReadBase + refBase;
        int qval = qualStr[readIdx];
        double qual = (useQual) ? sailfish::stringtools::phredToLogProb[qval] : sailfish::math::LOG_1;
        mismatchProfile[readIdx].increment(index, curReadBase, mass+p+qual);
        ++readIdx;
        ++transcriptIdx;
    }
    maxLen_ = std::max(maxLen_, static_cast<size_t>(read->core.l_qseq));
}

void ErrorModel::update(const ReadPair& hit, Transcript& ref, double p, double mass){
    if (mass == sailfish::math::LOG_0) { return; }
    bam1_t* leftRead = (hit.read1->core.mpos < hit.read2->core.mpos) ? hit.read1 : hit.read2;
    bam1_t* rightRead = (hit.read1->core.mpos < hit.read2->core.mpos) ? hit.read2 : hit.read1;

    if (leftRead) {
        update(leftRead, ref, p, mass, mismatchLeft_);
    }

    if (rightRead) {
        update(rightRead, ref, p, mass, mismatchRight_);
    }
}
