#ifndef __HDF5_WRITER_HPP__
#define __HDF5_WRITER_HPP__

#include <memory>
#include <mutex>

#include "spdlog/spdlog.h"

#include "SailfishSpinLock.hpp"
#include "SailfishOpts.hpp"
#include "ReadExperiment.hpp"
#include "EasyH5File.hpp"

class HDF5Writer {
  public:
    HDF5Writer(const boost::filesystem::path path, std::shared_ptr<spdlog::logger> logger);

    bool writeMeta(
	const SailfishOpts& opts, 
	const ReadExperiment& experiment);

    bool writeAbundances(
      const SailfishOpts& sopt,
      ReadExperiment& readExp);

    template <typename T>
    bool writeBootstrap(const std::vector<T>& estCounts);

   private:
     std::unique_ptr<EasyH5File> file_{nullptr};
     std::shared_ptr<spdlog::logger> logger_;
// only one writer thread at a time
#if defined __APPLE__
        spin_lock writeMutex_;
#else
        std::mutex writeMutex_;
#endif
        std::atomic<uint32_t> numBootstrapsWritten_{0};
};

#endif //__HDF5_WRITER_HPP__

