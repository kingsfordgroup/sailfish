#ifndef __SAILFISH_INDEX_VERSION_INFO_HPP__
#define __SAILFISH_INDEX_VERSION_INFO_HPP__

#include "spdlog/details/format.h"
#include "boost/filesystem.hpp"
#include "cereal/archives/json.hpp"

class SailfishIndexVersionInfo {
    public:
        /**
         * default constructor(s)
         */
        SailfishIndexVersionInfo() : indexVersion_(sailfish::indexVersion),
                                     kmerLength_(0) {}

        SailfishIndexVersionInfo(uint32_t indexVersionIn, uint32_t kmerLenIn):
            indexVersion_(indexVersionIn), kmerLength_(kmerLenIn) {}

        /**
         * Read the index version info from file
         */
        bool load(boost::filesystem::path& versionFile) {
            namespace bfs = boost::filesystem;
            if(!bfs::exists(versionFile)) {
                fmt::MemoryWriter infostr;
                infostr << "Error: The index version file " << versionFile.string()
                    << " doesn't seem to exist.  Please try re-building the sailfish "
                    "index.";
                throw std::invalid_argument(infostr.str());
            }
            std::ifstream ifs(versionFile.string());
            {
                cereal::JSONInputArchive iarchive(ifs); // Create an input archive
                iarchive(cereal::make_nvp("indexVersion", indexVersion_),
                        cereal::make_nvp("kmerLength", kmerLength_));
            }
            ifs.close();
            return true;
        }

        bool save(boost::filesystem::path& versionFile) {
            std::ofstream ofs(versionFile.string());
            {
                cereal::JSONOutputArchive oarchive(ofs);
                oarchive(cereal::make_nvp("indexVersion", indexVersion_),
                        cereal::make_nvp("kmerLength", kmerLength_));
            }
            ofs.close();
            return true;
        }

        uint32_t indexVersion() { return indexVersion_; }
        void indexVersion(uint32_t version) { indexVersion_ = version; }

        uint32_t kmerLength() { return kmerLength_; }
        void kmerLength(uint32_t len) { kmerLength_ = len; };

    private:
        uint32_t indexVersion_;
        uint32_t kmerLength_;
};

#endif // __SAILFISH_INDEX_VERSION_INFO_HPP__

