#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <array>
#include <unordered_map>
#include <limits>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "shark/Data/Dataset.h"
#include "shark/Algorithms/Trainers/PCA.h"
#include "shark/Algorithms/Trainers/RFTrainer.h"
#include "shark/Algorithms/Trainers/CARTTrainer.h"

#include "tensemble/TypeDef.h"
#include "tensemble/RandomForestRegressor.h"
#include "tensemble/RandomForestClassifier.h"
#include "tensemble/GBMRegressor.h"
#include "tensemble/GBMClassifier.h"
#include "tensemble/ReadData.h"

#include "CommonTypes.hpp"

#define DEFAULT_N_TREES 100
#define DEFAULT_N_JOBS 1
#define DEFAULT_MAX_FEATURES_RATIO 1.0
#define DEFAULT_MIN_SAMPLE_LEAF 5
#define DEFAULT_MAX_DEPTH   4
#define DEFAULT_SUBSAMPLE   1.0
#define DEFAULT_SPLIT_CRITERION CRITERION_MSE
#define DEFAULT_LOSS SQUARE_LOSS
#define DEFAULT_LEARN_RATE 0.1
#define DEFAULT_OOB 1
#define DEFAULT_VERBOSE 0
#define DEFAULT_BOOTSTRAP 1
#define DEFAULT_COMPUTE_IMPORTANCE 0

namespace bfs = boost::filesystem;
using Kmer = uint64_t;
using Sailfish::TranscriptFeatures;
using mpdec = boost::multiprecision::cpp_dec_float_100;

TranscriptFeatures parseFeature(std::ifstream& ifs) {
        TranscriptFeatures tf{};
        ifs >> tf.name;
        ifs >> tf.length;
        ifs >> tf.gcContent;
        for (auto i : boost::irange(size_t{0}, tf.diNucleotides.size())) {
                ifs >> tf.diNucleotides[i];
        }
        // eat the newline
        char junk;
        ifs.get(junk);
        return tf;
}

std::vector<TranscriptFeatures> parseFeatureFile(const bfs::path& featureFile) {
        std::ifstream ifile(featureFile.string());
        std::vector<TranscriptFeatures> feats;
        while (!ifile.eof()) {
                feats.emplace_back( parseFeature(ifile) );
                if (ifile.peek() == EOF) { break; }
        }
        ifile.close();
        return feats;
}

struct TranscriptResult{
        size_t length;
        double tpm;
        double rpkm;
        double kpkm;
        double approxCount;
};

struct ExpressionResults {
        std::vector<std::string> comments;
        std::unordered_map<std::string, TranscriptResult> expressions;
};

ExpressionResults parseSailfishFile(const bfs::path& expFile) {
        std::ifstream ifile(expFile.string());
        ExpressionResults res;
        while(!ifile.eof()) {

                if (ifile.peek() == '#') {
                        std::string comment;
                        std::getline(ifile, comment);
                        res.comments.emplace_back(comment);
                } else {
                        std::string tname;
                        TranscriptResult tr;
                        ifile >> tname;
                        ifile >> tr.length;
                        ifile >> tr.tpm;
                        ifile >> tr.rpkm;
                        ifile >> tr.kpkm;
                        ifile >> tr.approxCount;
                        res.expressions[tname] = tr;
                        // eat the newline
                        char nline; ifile.get(nline);
                }

                if (ifile.peek() == EOF) { break; }
        }

        return res;
}


void populateFromTPMs(vector<mpdec>& tpms,
                      vector<TranscriptFeatures>& features,
                      vector<size_t>& retainedRows,
                      ExpressionResults& sfres,
                      double estimatedReadLength,
                      double kmersPerRead,
                      uint64_t mappedKmers,
                      uint64_t merLen,
                      vector<mpdec>& rpkms,
                      vector<mpdec>& readCounts) {

  // compute the TPM normalization factor
  mpdec mpzero = 0;
  mpdec sumTPM = std::accumulate(tpms.begin(), tpms.end(), mpzero);
  mpdec norm = 1.0 / sumTPM;

  // compute the relative transcript fractions
  vector<mpdec> tfracs(tpms.size());
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    tfracs[i] = tpms[i] * norm;
  }

  // using the relative transcript fractions, compute the relative
  // nucleotide fractions (transcript fractions * length)
  vector<mpdec> tflens(tpms.size());
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    auto& name = features[i].name;
    auto& r = sfres.expressions[name];
    double l = (r.length - merLen + 1);
    tflens[i] = tfracs[i] * (l);
  }

  // normalize the nucleotide fractions
  mpdec tfnorm = 1.0 / std::accumulate(tflens.begin(), tflens.end(), mpzero);
  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    tflens[i] *= tfnorm;
  }

  uint64_t numReads = mappedKmers / kmersPerRead;

  // use the nucleotide fractions to compute the RPKMs
  double billion = pow(10,9);
  rpkms.clear(); rpkms.resize(features.size());
  readCounts.clear(); readCounts.resize(features.size());

  for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
    auto& name = features[i].name;
    auto& r = sfres.expressions[name];
    double l = (r.length - merLen + 1);
    rpkms[i] = billion * (tflens[i] / l);
    readCounts[i] = (tflens[i] * numReads) / kmersPerRead;
  }

}


int performBiasCorrection(
        bfs::path featureFile,
        bfs::path expressionFile,
        double estimatedReadLength,
        double kmersPerRead,
        uint64_t mappedKmers,
        uint32_t merLen,
        bfs::path outputFile,
        size_t numThreads) {
        //int argc, char* argv[]) {

        using shark::PCA;

        auto features = parseFeatureFile(featureFile);
        std::cerr << "parsed " << features.size() << " features\n";

        auto sfres = parseSailfishFile(expressionFile);
        std::cerr << "parsed " << sfres.expressions.size() << " expression values\n";

        std::vector<size_t> retainedRows;
        std::vector<double> retainedRPKMs;
        std::vector<std::string> retainedNames;

        double minLRPKM, maxLRPKM;
        minLRPKM = std::numeric_limits<double>::max();
        maxLRPKM = -minLRPKM;

        for (auto i : boost::irange(size_t{0}, features.size())) {
                auto& tname = features[i].name;
                auto rpkm = sfres.expressions[tname].kpkm; // ALTERATION
                shark::RealVector v(1);

                if ( rpkm >= 0.001 ) {
                        retainedRows.emplace_back(i);
                        retainedNames.push_back(tname);
                        v(0) = std::log(rpkm);
                        retainedRPKMs.push_back(v(0));
                        minLRPKM = std::min(minLRPKM, v(0));
                        maxLRPKM = std::max(maxLRPKM, v(0));
                }
        }


        //Pca pca;

        std::vector<shark::RealVector> featMat;
        std::vector<float> pcavec;
        size_t fnum = 0;
        for (auto r : retainedRows) {
                auto& f = features[r];
                shark::RealVector v(17);
                std::vector<double> fv(17);
                fv[0] = f.gcContent;
                v(0) = f.gcContent;
                pcavec.push_back(f.gcContent);
                for (auto i : boost::irange(size_t{0}, f.diNucleotides.size())) {
                        v(i+1) = f.diNucleotides[i];
                        fv[i+1] = f.diNucleotides[i];
                        pcavec.push_back(f.diNucleotides[i]);
                }
                featMat.emplace_back(v);
                ++fnum;
        }

        /*
        std::cerr << "solving . . . ";
        pca.Calculate(pcavec, retainedRows.size(), 17, true, false, false);
        std::cerr << "done\n";
        size_t numRetained = pca.thresh95();
        std::cerr << "numRetained = " << numRetained << "\n";
        std::vector<float> sdevs = pca.sd();
        std::cerr << "sdevs.size() " << sdevs.size() << "\n";
        */

        shark::UnlabeledData<shark::RealVector> Xsub = shark::createDataFromRange(featMat);

        PCA pcao(Xsub, true);
        auto evals = pcao.eigenvalues();
        double totalVariance = 0.0;
        for ( auto e : evals ) { std::cerr << e << "\n"; totalVariance += e; }
        std::cerr << "totalVariance: " << totalVariance << "\n";
        double varCutoff = 0.95;
        double varSum = 0.0;
        size_t dimCutoff = 0;
        size_t currentDim = 0;
        for ( auto e : evals ) {
                ++currentDim;
                varSum += e;
                std::cerr << "ev: " << e <<  "\n";
                if (varSum / totalVariance >= varCutoff) {
                        dimCutoff = currentDim;
                        break;
                }

        }
        std::cerr << varCutoff * 100.0 << "% of the variance is explained by " << dimCutoff << " dimensions\n";

        /*
        std::vector<float> scores = pca.scores();
        std::cerr << "score.size() " << scores.size() << "\n";
        std::cerr << "pts * dimes = " << retainedRows.size() * dimCutoff << "\n";
        */

        shark::LinearModel<> enc;
        pcao.encoder(enc, dimCutoff);
        auto encodedXsub = enc(Xsub);

        //shark::UnlabeledData<shark::RealVector> X(18, features.size());
        Data train;
        train.set_size(retainedRows.size(), dimCutoff+1);

        //size_t le = 0;
        size_t c = 0;
        for (auto r : retainedRows) {
                train.X[c][0] = std::log(static_cast<double>(features[r].length));

                for (auto j : boost::irange(size_t{1}, dimCutoff+1)) {
                        train.X[c][j] = encodedXsub.element(c)(j-1);
                        //train.X[c][j] = scores[le];
                        //train2.X[c][j] = scores[le];
                        //++le;
                        //std::cerr << "Train [" << c <<"][" << j        << "] = " << train.X[c][j] << "\n";
                }
                //le += (17 - dimCutoff);
                train.y[c] = retainedRPKMs[c];// std::log(sfres[features[r].name].expressions.rpkm);
                ++c;
        }

        /** Random Forest Regression **/

        auto reg = std::unique_ptr<RandomForestRegressor>(new RandomForestRegressor(
                500,
                train.n_features,
                4, // max tree depth
                1, // min_samples_leaf
                1.0, // features ratio
                true, // bootstrap
                true, //out-of-bag
                true, // compute importance
                0, // random seed
                numThreads, // num jobs
                true // verbose
        ));


        /*
        auto reg = std::unique_ptr<GBMRegressor>(new GBMRegressor(
                SQUARE_LOSS,
                500,
                train.n_features,
                4, // max tree depth
                1, //min sample leaf
                1.0, // max features ratio
                0.5, // subsample
                0.05, //learn rate
                true, //out-of-bag
                true, // compute imporance
                34239, // random seed
                numThreads, // num jobs
                true // verbose
        ));
        */

        std::cerr << "there are " << train.n_samples << " samples\n";
        std::cerr << "there are " << train.n_features << " features\n";
        reg->build(train.X, train.y, train.n_samples);

        std::vector<REAL> pred(train.n_samples, 0.0);
        reg->predict(train.X, &pred[0], train.n_samples, train.n_features);

  REAL trn_rmse=rmse(&pred[0], train.y, train.n_samples);
  REAL trn_r2=R2(&pred[0], train.y, train.n_samples);
  std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";

        double grandMean = 0.0;
        for (auto i : boost::irange(size_t{0}, size_t{train.n_samples})) {
                grandMean += retainedRPKMs[i];
        }
        grandMean /= train.n_samples;

        for (auto i : boost::irange(size_t{0}, size_t{train.n_samples})) {
                pred[i] = grandMean + (retainedRPKMs[i] - pred[i]);
        }

        /**
        for (auto i : boost::irange(size_t{0}, size_t{train.n_samples})) {
                pred[i] = retainedRPKMs[i] - pred[i];
        }

        auto mmpred = std::minmax_element(pred.begin(), pred.end());
        double minPred = *(mmpred.first);
        double maxPred = *(mmpred.second);

        double scale = std::fabs(maxLRPKM - minLRPKM) / std::fabs(maxPred - minPred);

        std::cerr << "min,max LRPKM : " << minLRPKM << ", " << maxLRPKM << "\n";
        std::cerr << "min, max pred : " << minPred << ", " << maxPred  << "\n";
        std::cerr << "SCALE: " << scale << "\n";


        minPred = std::numeric_limits<double>::max();
        for (auto i : boost::irange(size_t{0}, size_t{train.n_samples})) {
                pred[i] *= scale;
                minPred = std::min(minPred, pred[i]);
        }

        double shift{minLRPKM - minPred};
        minPred = std::numeric_limits<double>::max();
        maxPred = -minPred;
        for (auto i : boost::irange(size_t{0}, size_t{train.n_samples})) {
                pred[i] += shift;
                minPred = std::min(minPred, pred[i]);
                maxPred = std::max(maxPred, pred[i]);
        }
        **/
   trn_rmse=rmse(&pred[0], train.y, train.n_samples);
   trn_r2=R2(&pred[0], train.y, train.n_samples);
  std::cerr << "Train RMSE=" << trn_rmse << ", Correlation Coefficient=" << trn_r2 << "\n";


        //shark::UnlabeledData<shark::RealVector> X()

        std::ofstream ofile(outputFile.string());
        for (auto& c : sfres.comments) {
                ofile << c << "\n";
        }


        size_t retainedCnt = 0;
        vector<mpdec> kpkms(features.size());
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = sfres.expressions[name];
          if (i == retainedRows[retainedCnt]) {
            kpkms[i] = std::exp(pred[retainedCnt]);
            ++retainedCnt;
          } else {
            kpkms[i] = r.kpkm;
          }
        }


        // compute estimated TPM from the KPKMS
        mpdec mpzero = 0;
        // normalize the KPKMS --- these will estimate the tau_i
        mpdec sumKPKM = std::accumulate(kpkms.begin(), kpkms.end(), mpzero);
        mpdec norm = 1.0 / sumKPKM;

        // then multiply by 10^6 to get TPM_i
        double million = pow(10, 6);
        vector<mpdec> tpms(kpkms.size());
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          tpms[i] = kpkms[i] * norm * million;
        }

        /*
        size_t retainedCnt = 0;
        vector<mpdec> tpms(features.size());
        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = sfres.expressions[name];
          if (i == retainedRows[retainedCnt]) {
            tpms[i] = std::exp(pred[retainedCnt]);
            ++retainedCnt;
          } else {
            tpms[i] = r.tpm;
          }
        }
        */

        vector<mpdec> rpkms;
        vector<mpdec> readCounts;
        populateFromTPMs(tpms, features, retainedRows,
                         sfres, estimatedReadLength, kmersPerRead,
                         mappedKmers, merLen, rpkms, readCounts);

        for (auto i : boost::irange(size_t{0}, size_t{features.size()})) {
          auto& name = features[i].name;
          auto& r = sfres.expressions[name];
          auto length = r.length;
          ofile << name << '\t' << r.length << '\t' << tpms[i] << '\t'
                << ((length - estimatedReadLength + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                << ((length - merLen + 1) > 0 ? kpkms[i] : 0.0) << '\t'
                << readCounts[i] << '\n';
        }
        std::cerr << "retainedCnt = " << retainedCnt << ", nsamps = " << train.n_samples << "\n";

        ofile.close();

}
