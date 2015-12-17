#include <ctime>
#include "HDF5Writer.hpp"

#include "EasyH5File.hpp"
#include "SailfishOpts.hpp"

HDF5Writer::HDF5Writer(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger) :
	               logger_(logger) {
  file_.reset(new EasyH5File(path.string()));
}

bool HDF5Writer::writeMeta(
    const SailfishOpts& opts, 
    const ReadExperiment& experiment) {

  file_->addGroup("/aux");
  auto numBootstraps = opts.numBootstraps;
  auto numSamples = (numBootstraps > 0) ? numBootstraps : opts.numGibbsSamples;
  if (numSamples > 0) {
    file_->addGroup("/bootstrap");
  }

  std::string sampType = "bootstrap";
  if (numBootstraps == 0 and numSamples > 0) {
    sampType = "gibbs";
  }
  std::vector<std::string> sampTypeVec{sampType};
  file_->writeVectorToGroup(sampTypeVec, "/aux", "sample_type");

  std::vector<uint64_t> numSampVec{experiment.numObservedFragments()};
  file_->writeVectorToGroup(numSampVec, "/aux", "num_processed");

  // For spoofing kallisto version
  std::vector<std::string> kalVer{"0.42.4"}; 
  file_->writeVectorToGroup(kalVer, "/aux", "kallisto_version");


  std::vector<int> kalIndexVersion{10};
  file_->writeVectorToGroup(kalIndexVersion, "/aux", "index_version");

  std::vector<std::string> call{"quant"};
  file_->writeVectorToGroup(call, "/aux", "call");


  std::time_t result = std::time(NULL);
  auto tstring = std::asctime(std::localtime(&result));
  std::vector<std::string> startTime{tstring};
  file_->writeVectorToGroup(call, "/aux", "start_time");
  return true;
}

template <typename T>
bool HDF5Writer::writeBootstrap(const std::vector<T>& abund) {
#if defined __APPLE__
            spin_lock::scoped_lock sl(writeMutex_);
#else
            std::lock_guard<std::mutex> lock(writeMutex_);
#endif
	    std::string bsStr = "bs" + std::to_string(numBootstrapsWritten_.load());
	    file_->writeVectorToGroup(abund, "/bootstrap", bsStr);
            logger_->info("wrote {} bootstraps", numBootstrapsWritten_.load()+1);
            ++numBootstrapsWritten_;
            return true;
}

bool HDF5Writer::writeAbundances(
    const SailfishOpts& sopt,
    ReadExperiment& readExp) {

  bool useScaledCounts = (sopt.allowOrphans == false);
  double numMappedFrags = readExp.numMappedFragments();

  std::vector<Transcript>& transcripts_ = readExp.transcripts();
  for (auto& transcript : transcripts_) {
    transcript.projectedCounts = useScaledCounts ?
      (transcript.mass() * numMappedFrags) : transcript.estCount();
  }

  double tfracDenom{0.0};
  for (auto& transcript : transcripts_) {
    double refLength = sopt.noEffectiveLengthCorrection ?
      transcript.RefLength :
      transcript.EffectiveLength;
      tfracDenom += (transcript.projectedCounts / numMappedFrags) / refLength;
  }

  double million = 1000000.0;

  std::vector<double> estCounts;
  estCounts.reserve(transcripts_.size());
  
  std::vector<std::string> names;
  names.reserve(transcripts_.size());
 
  std::vector<double> effectiveLengths;
  effectiveLengths.reserve(transcripts_.size());
 
  std::vector<uint32_t> lengths;
  lengths.reserve(transcripts_.size());

  std::vector<double> tpms;
  tpms.reserve(transcripts_.size());
 
  for (auto& txp : transcripts_) {
    double count = txp.projectedCounts;
    double effLen = sopt.noEffectiveLengthCorrection ? 
	txp.RefLength : txp.EffectiveLength;
    double npm = (count / numMappedFrags);
    double tfrac = (npm / effLen) / tfracDenom;
    double tpm = tfrac * million;
    
    estCounts.push_back(count);
    lengths.push_back(txp.RefLength);
    effectiveLengths.push_back(txp.EffectiveLength);
    names.push_back(txp.RefName);


    tpms.push_back(tpm);
  }

  file_->writeVectorToGroup(estCounts, "/", "est_counts");
  file_->writeVectorToGroup(names, "/aux", "est_counts");
  file_->writeVectorToGroup(effectiveLengths, "/aux", "eff_lengths");
  file_->writeVectorToGroup(lengths, "/aux", "lengths");
  file_->writeVectorToGroup(tpms, "/aux", "tpms");
  return true;
}

template 
bool HDF5Writer::writeBootstrap<double>(const std::vector<double>& abund);

template 
bool HDF5Writer::writeBootstrap<int>(const std::vector<int>& abund);

