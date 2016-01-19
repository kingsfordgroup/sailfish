#ifndef __GZIP_WRITER_HPP__
#define __GZIP_WRITER_HPP__

#include <memory>
#include <mutex>

#include "spdlog/spdlog.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "SailfishSpinLock.hpp"
#include "SailfishOpts.hpp"
#include "ReadExperiment.hpp"

class GZipWriter {
  public:
    GZipWriter(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger);

    ~GZipWriter();

    bool writeEquivCounts(
        const SailfishOpts& opts,
        ReadExperiment& experiment);

    bool writeMeta(
	const SailfishOpts& opts,
	const ReadExperiment& experiment,
        const std::string& tstring  = "now"  // the start time of the run
	);

    bool writeAbundances(
      const SailfishOpts& sopt,
      ReadExperiment& readExp);

    template <typename T>
    bool writeBootstrap(const std::vector<T>& abund);

   private:
     boost::filesystem::path path_;
     boost::filesystem::path bsPath_;
     std::shared_ptr<spdlog::logger> logger_;
     std::unique_ptr<boost::iostreams::filtering_ostream> bsStream_{nullptr};
// only one writer thread at a time
#if defined __APPLE__
        spin_lock writeMutex_;
#else
        std::mutex writeMutex_;
#endif
        std::atomic<uint32_t> numBootstrapsWritten_{0};
};

#endif //__GZIP_WRITER_HPP__
