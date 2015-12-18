#ifndef __SAILFISH_INDEX_HPP__
#define __SAILFISH_INDEX_HPP__

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "spdlog/spdlog.h"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"

#include "RapMapSAIndex.hpp"
#include "IndexHeader.hpp"
#include "SailfishConfig.hpp"
#include "SailfishIndexVersionInfo.hpp"

// declaration of quasi index function
int rapMapSAIndex(int argc, char* argv[]);

class SailfishIndex{
    public:
        SailfishIndex(std::shared_ptr<spdlog::logger>& logger) :
            loaded_(false), versionInfo_(sailfish::indexVersion, 0),
            logger_(logger) {}


        void load(const boost::filesystem::path& indexDir) {
            namespace bfs = boost::filesystem;

            // Check if version file exists and, if so, read it.
            boost::filesystem::path versionPath = indexDir / "versionInfo.json";
            versionInfo_.load(versionPath);
            if (versionInfo_.indexVersion() < sailfish::indexVersion) {
                fmt::MemoryWriter infostr;
                infostr << "[Error]: The index version appears to be too old. "
                        << "Please re-build the sailfish index.";
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
        bool is64BitQuasi() { return largeIndex_; }
        RapMapSAIndex<int32_t>* quasiIndex32() { return quasiIndex32_.get(); }
        RapMapSAIndex<int64_t>* quasiIndex64() { return quasiIndex64_.get(); }

        const char* transcriptomeSeq() {
            if (loaded_) {
                if (is64BitQuasi()) {
                    return quasiIndex64_->seq.c_str();
                } else {
                    return quasiIndex32_->seq.c_str();
                }
            } else {
                return nullptr;
            }
        }

        uint64_t transcriptOffset(uint64_t id) {
            if (loaded_) {
                if (is64BitQuasi()) {
                    return quasiIndex64_->txpOffsets[id];
                } else {
                    return quasiIndex32_->txpOffsets[id];
                }
            } else {
                return std::numeric_limits<uint64_t>::max();
            }
        }

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

                // Read the index header and determine
                // the size.
                IndexHeader h;
                std::ifstream indexStream(indexStr + "header.json");
                {
                    cereal::JSONInputArchive ar(indexStream);
                    ar(h);
                }
                indexStream.close();

                if (h.bigSA()) {
                    largeIndex_ = true;
                    fmt::print(stderr, "Loading 64-bit quasi index");
                    quasiIndex64_.reset(new RapMapSAIndex<int64_t>);
                    if (!quasiIndex64_->load(indexStr)) {
                        fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                        fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                        std::exit(1);
                    }
                } else {
                    fmt::print(stderr, "Loading 32-bit quasi index");
                    quasiIndex32_.reset(new RapMapSAIndex<int32_t>);
                    if(!quasiIndex32_->load(indexStr)) {
                        fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                        fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                        std::exit(1);
                    }
                }
            }
            logger_->info("done");
            return true;
        }


        bool loaded_;
        SailfishIndexVersionInfo versionInfo_;
        //std::unique_ptr<RapMapSAIndex> quasiIndex_{nullptr};
        // Can't think of a generally better way to do this now
        // without making the entire code-base look crazy
        bool largeIndex_{false};
        std::unique_ptr<RapMapSAIndex<int32_t>> quasiIndex32_{nullptr};
        std::unique_ptr<RapMapSAIndex<int64_t>> quasiIndex64_{nullptr};
        std::shared_ptr<spdlog::logger> logger_;
};

#endif // __SAILFISH_INDEX_HPP__
