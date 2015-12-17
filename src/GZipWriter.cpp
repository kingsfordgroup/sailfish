#include <ctime>
#include <fstream>

#include "cereal/archives/json.hpp"

#include "GZipWriter.hpp"
#include "SailfishOpts.hpp"

GZipWriter::GZipWriter(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger) :
  path_(path), logger_(logger) {
}

GZipWriter::~GZipWriter() {
  if (bsStream_) {
    bsStream_->reset();
  }
}

template <typename T>
bool writeVectorToFile(boost::filesystem::path path,
                       const std::vector<T>& vec) {

    {
        bool binary = std::is_same<T, std::string>::value;
        auto flags = std::ios_base::out | std::ios_base::binary;

        boost::iostreams::filtering_ostream out;
        out.push(boost::iostreams::gzip_compressor(6));
        out.push(boost::iostreams::file_sink(path.string(), flags));

        size_t num = vec.size();
        size_t elemSize = sizeof(typename std::vector<T>::value_type);
        // We have to get rid of constness below, but this should be OK
        out.write(reinterpret_cast<char*>(const_cast<T*>(vec.data())),
                  num * elemSize);
        out.reset();
    }
    return true;
}


bool GZipWriter::writeMeta(
    const SailfishOpts& opts,
    const ReadExperiment& experiment) {

  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / "aux";
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  auto numBootstraps = opts.numBootstraps;
  auto numSamples = (numBootstraps > 0) ? numBootstraps : opts.numGibbsSamples;
  if (numSamples > 0) {
      bsPath_ = auxDir / "bootstrap";
      //bfs::path bsDir = auxDir / "bootstrap";
      bool bsSuccess = boost::filesystem::create_directories(bsPath_);
      {

          boost::iostreams::filtering_ostream nameOut;
          nameOut.push(boost::iostreams::gzip_compressor(6));
          auto bsFilename = bsPath_ / "names.tsv.gz";
          nameOut.push(boost::iostreams::file_sink(bsFilename.string(), std::ios_base::out));

          auto& transcripts = experiment.transcripts();
          size_t numTxps = transcripts.size();
          if (numTxps == 0) { return false; }
          for (size_t tn = 0; tn < numTxps; ++tn) {
              auto& t  = transcripts[tn];
              nameOut << t.RefName;
              if (tn < numTxps - 1) {
                  nameOut << '\t';
              }
          }
          nameOut << '\n';
          nameOut.reset();
      }

  }

  bfs::path fldPath = auxDir / "fld.gz";
  writeVectorToFile(fldPath, experiment.fragLengthDist());

  bfs::path normBiasPath = auxDir / "expected_bias.gz";
  writeVectorToFile(normBiasPath, experiment.expectedBias());

  bfs::path obsBiasPath = auxDir / "observed_bias.gz";
  const auto& bcounts = experiment.readBias().counts;
  std::vector<int32_t> observedBias(bcounts.size(), 0);
  std::copy(bcounts.begin(), bcounts.end(), observedBias.begin());
  writeVectorToFile(obsBiasPath, observedBias);

  bfs::path info = auxDir / "meta_info.json";

  {
      std::ofstream os(info.string());
      cereal::JSONOutputArchive oa(os);

      std::string sampType = "none";
      if (numBootstraps == 0 and numSamples > 0) {
          sampType = "gibbs";
      }
      if (numBootstraps > 0) {
          sampType = "bootstrap";
      }

      auto& transcripts = experiment.transcripts();
      oa(cereal::make_nvp("sf_version", std::string(sailfish::version)));
      oa(cereal::make_nvp("samp_type", sampType));
      oa(cereal::make_nvp("frag_dist_length", experiment.fragLengthDist().size()));
      oa(cereal::make_nvp("num_bias_bins", bcounts.size()));
      oa(cereal::make_nvp("num_targets", transcripts.size()));
      oa(cereal::make_nvp("num_bootstraps", numBootstraps));
      oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
      oa(cereal::make_nvp("call", std::string("quant")));

      std::time_t result = std::time(NULL);
      std::string tstring(std::asctime(std::localtime(&result)));

      oa(cereal::make_nvp("start_time", tstring));
  }
  // For spoofing kallisto version
  //std::vector<std::string> kalVer{"0.42.4"};
  //file_->writeVectorToGroup(kalVer, "/aux", "kallisto_version");

  //std::vector<int> kalIndexVersion{10};
  //file_->writeVectorToGroup(kalIndexVersion, "/aux", "index_version");
  return true;
}

bool GZipWriter::writeAbundances(
    const SailfishOpts& sopt,
    ReadExperiment& readExp) {
  return false;
  /*
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
  */
}

template <typename T>
bool GZipWriter::writeBootstrap(const std::vector<T>& abund) {
#if defined __APPLE__
            spin_lock::scoped_lock sl(writeMutex_);
#else
            std::lock_guard<std::mutex> lock(writeMutex_);
#endif
	    if (!bsStream_) {
	      bsStream_.reset(new boost::iostreams::filtering_ostream);
	      bsStream_->push(boost::iostreams::gzip_compressor(6));
	      auto bsFilename = bsPath_ / "bootstraps.gz";
	      bsStream_->push(
                  boost::iostreams::file_sink(bsFilename.string(),
                                              std::ios_base::out | std::ios_base::binary));
	    }

	    boost::iostreams::filtering_ostream& ofile = *bsStream_;
	    size_t num = abund.size();
        size_t elSize = sizeof(typename std::vector<T>::value_type);
        ofile.write(reinterpret_cast<char*>(const_cast<T*>(abund.data())),
                    elSize * num);
        /*
        for (size_t tn = 0; tn < num; ++tn) {
            auto& a  = abund[tn];
            ofile << a;
            if (tn < num - 1) {
                ofile << '\t';
            }
        }
        ofile << '\n';
        */
        logger_->info("wrote {} bootstraps", numBootstrapsWritten_.load()+1);
        ++numBootstrapsWritten_;
        return true;
}

template
bool GZipWriter::writeBootstrap<double>(const std::vector<double>& abund);

template
bool GZipWriter::writeBootstrap<int>(const std::vector<int>& abund);

