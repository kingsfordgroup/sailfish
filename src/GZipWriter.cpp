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

/**
 * Creates a new gzipped file (path) and writes the contents
 * of the vector (vec) to the file in binary.
 */
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

/**
 * Write the ``main'' metadata to file.  Currently this includes:
 *   -- Names of the target id's if bootstrapping / gibbs is performed
 *   -- The fragment length distribution
 *   -- The expected and observed bias values
 *   -- A json file with information about the run
 */
bool GZipWriter::writeMeta(
    const SailfishOpts& opts,
    const ReadExperiment& experiment,
    const std::string& tstring // the start time of the run
    ) {

  namespace bfs = boost::filesystem;

  bfs::path auxDir = path_ / "aux";
  bool auxSuccess = boost::filesystem::create_directories(auxDir);

  auto numBootstraps = opts.numBootstraps;
  auto numSamples = (numBootstraps > 0) ? numBootstraps : opts.numGibbsSamples;
  if (numSamples > 0) {
      bsPath_ = auxDir / "bootstrap";
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
      oa(cereal::make_nvp("bias_correct", opts.biasCorrect));
      oa(cereal::make_nvp("num_bias_bins", bcounts.size()));
      oa(cereal::make_nvp("num_targets", transcripts.size()));
      oa(cereal::make_nvp("num_bootstraps", numBootstraps));
      oa(cereal::make_nvp("num_processed", experiment.numObservedFragments()));
      oa(cereal::make_nvp("num_mapped", experiment.numMappedFragments()));
      oa(cereal::make_nvp("percent_mapped", experiment.mappingRate() * 100.0));
      oa(cereal::make_nvp("call", std::string("quant")));
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
  
  using sailfish::math::LOG_0;
  using sailfish::math::LOG_1;
  namespace bfs = boost::filesystem;

  bool useScaledCounts = (sopt.allowOrphans == false);
  bfs::path fname = path_ / "quant.sf";

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

  /* No comments for now
  if (headerComments.length() > 0) {
    fmt::print(output.get(), "{}", headerComments);
  }
  */

  // The header
  fmt::print(output.get(), "Name\tLength\tEffectiveLength\tTPM\tNumReads\n");

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
  // Now posterior has the transcript fraction
  for (auto& transcript : transcripts_) {
    auto effLen = sopt.noEffectiveLengthCorrection ?
      transcript.RefLength :
      transcript.EffectiveLength;
    double count = transcript.projectedCounts;
    double npm = (transcript.projectedCounts / numMappedFrags);
    double tfrac = (npm / effLen) / tfracDenom;
    double tpm = tfrac * million;
    fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
	transcript.RefName, transcript.RefLength, effLen,
	tpm, count);
  }

  return true;
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

