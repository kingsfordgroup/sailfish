#ifndef __SAILFISH_INDEX_HPP__
#define __SAILFISH_INDEX_HPP__

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "spdlog/spdlog.h"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"

#include "RapMapSAIndex.hpp"
#include "SailfishConfig.hpp"
#include "SailfishIndexVersionInfo.hpp"

// declaration of quasi index function
int rapMapSAIndex(int argc, char* argv[]);

class SailfishIndex{
    public:
        SailfishIndex(std::shared_ptr<spdlog::logger>& logger) :
            loaded_(false), versionInfo_(0, 0), logger_(logger) {}


        void load(const boost::filesystem::path& indexDir) {
            namespace bfs = boost::filesystem;

            // Check if version file exists and, if so, read it.
            boost::filesystem::path versionPath = indexDir / "versionInfo.json";
            versionInfo_.load(versionPath);
            if (versionInfo_.indexVersion() == 0) {
                fmt::MemoryWriter infostr;
                infostr << "Error: The index version file " << versionPath.string()
                    << " doesn't seem to exist.  Please try re-building the salmon "
                    "index.";
                throw std::invalid_argument(infostr.str());
            }
            // Check index version compatibility here
            loadQuasiIndex_(indexDir);
            loaded_ = true;
        }

        bool build(boost::filesystem::path indexDir,
                std::vector<const char*>& argVec, uint32_t k) {
            return buildQuasiIndex_(indexDir, argVec, k);
        }

        bool loaded() { return loaded_; }
        RapMapSAIndex* quasiIndex() { return quasiIndex_.get(); }

    private:
        bool buildQuasiIndex_(boost::filesystem::path indexDir,
                             std::vector<const char*>& quasiArgVec,
                             uint32_t k) {
            namespace bfs = boost::filesystem;
            char* quasiArgv[] = { const_cast<char*>(quasiArgVec[0]),
                const_cast<char*>(quasiArgVec[1]),
                const_cast<char*>(quasiArgVec[2]),
                const_cast<char*>(quasiArgVec[3]),
                const_cast<char*>(quasiArgVec[4]),
                const_cast<char*>(quasiArgVec[5]),
                const_cast<char*>(quasiArgVec[6])
            };
            int quasiArgc = 7;

            int ret = rapMapSAIndex(quasiArgc, quasiArgv);

            bfs::path versionFile = indexDir / "versionInfo.json";
            versionInfo_.indexVersion(sailfish::indexVersion);
            versionInfo_.kmerLength(k);
            versionInfo_.save(versionFile);
            return (ret == 0);
        }

        bool loadQuasiIndex_(const boost::filesystem::path& indexDir) {
            namespace bfs = boost::filesystem;
            logger_->info("Loading Quasi index");
            // Read the actual Quasi index
            { // quasi-based
                boost::filesystem::path indexPath = indexDir;
                std::string indexStr = indexDir.string();
                if (indexStr.back() != '/') { indexStr.push_back('/'); }
                quasiIndex_.reset(new RapMapSAIndex);
                if (!quasiIndex_->load(indexStr)) {
                    fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                    fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                    std::exit(1);
                }
            }
            logger_->info("done");
            return true;
        }


        bool loaded_;
        SailfishIndexVersionInfo versionInfo_;
        std::unique_ptr<RapMapSAIndex> quasiIndex_{nullptr};
        std::shared_ptr<spdlog::logger> logger_;
};

#endif // __SAILFISH_INDEX_HPP__
