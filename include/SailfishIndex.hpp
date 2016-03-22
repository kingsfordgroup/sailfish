#ifndef __SAILFISH_INDEX_HPP__
#define __SAILFISH_INDEX_HPP__

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "spdlog/spdlog.h"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"

#include "MPHMap.hpp"
#include "RapMapSAIndex.hpp"
#include "IndexHeader.hpp"
#include "SailfishConfig.hpp"
#include "SailfishIndexVersionInfo.hpp"

// declaration of quasi index function
int rapMapSAIndex(int argc, char* argv[]);

template <typename IndexT> 
using DenseHash = google::dense_hash_map<uint64_t, 
                                         rapmap::utils::SAInterval<IndexT>, 
                                         rapmap::utils::KmerKeyHasher>;
template <typename IndexT> 
using PerfectHash = MPHMap<uint64_t, 
                           std::pair<uint64_t, 
                                     rapmap::utils::SAInterval<IndexT>>>;

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
		   std::vector<std::string>& argVec, uint32_t k) {
            return buildQuasiIndex_(indexDir, argVec, k);
        }

        bool loaded() { return loaded_; }
        bool is64BitQuasi() { return largeIndex_; }
        bool isPerfectHashQuasi() { return perfectHashIndex_; }

        RapMapSAIndex<int32_t, DenseHash<int32_t>>* quasiIndex32() { return quasiIndex32_.get(); }
        RapMapSAIndex<int64_t, DenseHash<int64_t>>* quasiIndex64() { return quasiIndex64_.get(); }

        RapMapSAIndex<int32_t, PerfectHash<int32_t>>* quasiIndexPerfectHash32() { return quasiIndexPerfectHash32_.get(); }
        RapMapSAIndex<int64_t, PerfectHash<int64_t>>* quasiIndexPerfectHash64() { return quasiIndexPerfectHash64_.get(); }


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
			      std::vector<std::string>& quasiArgVec,
                             uint32_t k) {
            namespace bfs = boost::filesystem;
	    int quasiArgc = static_cast<int>(quasiArgVec.size());
	    char** quasiArgv = new char*[quasiArgc];
	    for (size_t i = 0; i < quasiArgc; ++i) {
	      auto& arg = quasiArgVec[i];
	      quasiArgv[i] = new char[arg.size() + 1];
	      std::strcpy(quasiArgv[i], arg.c_str());
	    }
	    /*
            char* quasiArgv[] = { const_cast<char*>(quasiArgVec[0]),
                const_cast<char*>(quasiArgVec[1]),
                const_cast<char*>(quasiArgVec[2]),
                const_cast<char*>(quasiArgVec[3]),
                const_cast<char*>(quasiArgVec[4]),
                const_cast<char*>(quasiArgVec[5]),
                const_cast<char*>(quasiArgVec[6])
            };
            int quasiArgc = 7;
	    */

            int ret = rapMapSAIndex(quasiArgc, quasiArgv);

            bfs::path versionFile = indexDir / "versionInfo.json";
            versionInfo_.indexVersion(sailfish::indexVersion);
            versionInfo_.kmerLength(k);
            versionInfo_.save(versionFile);

	    // Free the memory used for the arg vector
	    for (size_t i = 0; i < quasiArgc; ++i) {
	      delete quasiArgv[i];
	    }
	    delete [] quasiArgv;
	    
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

		// Is the index using a perfect hash
		perfectHashIndex_ = h.perfectHash();
		
                if (h.bigSA()) {
                    largeIndex_ = true;
                    fmt::print(stderr, "Loading 64-bit quasi index");
		    if (perfectHashIndex_) {
		      quasiIndexPerfectHash64_.reset(new RapMapSAIndex<int64_t, PerfectHash<int64_t>>);
		      if (!quasiIndexPerfectHash64_->load(indexStr)) {
			fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
			fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
			std::exit(1);
		      }
		    } else {
		      quasiIndex64_.reset(new RapMapSAIndex<int64_t, DenseHash<int64_t>>);
		      if (!quasiIndex64_->load(indexStr)) {
			fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
			fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
			std::exit(1);
		      }
		    }
		} else { // 32-bit index
                    fmt::print(stderr, "Loading 32-bit quasi index");
		    
		    if (perfectHashIndex_) {
		      quasiIndexPerfectHash32_.reset(new RapMapSAIndex<int32_t, PerfectHash<int32_t>>);
		      if (!quasiIndexPerfectHash32_->load(indexStr)) {
			fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
			fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
			std::exit(1);
		      }
		    } else {
		      quasiIndex32_.reset(new RapMapSAIndex<int32_t, DenseHash<int32_t>>);
		      if (!quasiIndex32_->load(indexStr)) {
			fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
			fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
			std::exit(1);
		      }
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
        bool perfectHashIndex_{false};
        std::unique_ptr<RapMapSAIndex<int32_t, DenseHash<int32_t>>> quasiIndex32_{nullptr};
        std::unique_ptr<RapMapSAIndex<int64_t, DenseHash<int64_t>>> quasiIndex64_{nullptr};
        std::unique_ptr<RapMapSAIndex<int32_t, PerfectHash<int32_t>>> quasiIndexPerfectHash32_{nullptr};
        std::unique_ptr<RapMapSAIndex<int64_t, PerfectHash<int64_t>>> quasiIndexPerfectHash64_{nullptr};
        std::shared_ptr<spdlog::logger> logger_;
};

#endif // __SAILFISH_INDEX_HPP__
