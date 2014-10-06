#include "SalmonUtils.hpp"
#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include "SailfishMath.hpp"
#include "LibraryFormat.hpp"
#include "ReadExperiment.hpp"

namespace salmon {
namespace utils {

    template <typename ExpLib>
    void writeAbundances(ExpLib& alnLib,
                         boost::filesystem::path& fname,
                         std::string headerComments) {
        using sailfish::math::LOG_0;
        using sailfish::math::LOG_1;

        std::unique_ptr<std::FILE, int (*)(std::FILE *)> output(std::fopen(fname.c_str(), "w"), std::fclose);

        fmt::print(output.get(), "{}", headerComments);
        fmt::print(output.get(), "# Name\tLength\tTPM\tFPKM\tNumReads\n");


        auto& refs = alnLib.transcripts();
        auto numMappedReads = alnLib.numMappedReads();
        const double logBillion = std::log(1000000000.0);
        const double million = 1000000.0;
        const double logNumFragments = std::log(static_cast<double>(numMappedReads));
        auto clusters = alnLib.clusterForest().getClusters();
        size_t clusterID = 0;
        for(auto cptr : clusters) {

            double logClusterMass = cptr->logMass();
            double logClusterCount = std::log(static_cast<double>(cptr->numHits()));

            if (logClusterMass == LOG_0) {
                std::cerr << "Warning: cluster " << clusterID << " has 0 mass!\n";
            }

            bool requiresProjection{false};

            auto& members = cptr->members();
            size_t clusterSize{0};
            for (auto transcriptID : members) {
                Transcript& t = refs[transcriptID];
                t.uniqueCounts = t.uniqueCount();
                t.totalCounts = t.totalCount();
                //clusterCount += t.totalCounts;
            }

            for (auto transcriptID : members) {
                Transcript& t = refs[transcriptID];
                double logTranscriptMass = t.mass(false);
                if (logTranscriptMass == LOG_0) {
                    t.projectedCounts = 0;
                } else {
                    double logClusterFraction = logTranscriptMass - logClusterMass;
                    t.projectedCounts = std::exp(logClusterFraction + logClusterCount);
                    requiresProjection |= t.projectedCounts > static_cast<double>(t.totalCounts) or
                        t.projectedCounts < static_cast<double>(t.uniqueCounts);
                }
                ++clusterSize;
            }

            if (clusterSize > 1 and requiresProjection) {
                cptr->projectToPolytope(refs);
            }
            ++clusterID;
        }

        auto& transcripts_ = refs;
        double tfracDenom{0.0};
        for (auto& transcript : transcripts_) {
            tfracDenom += (transcript.projectedCounts / numMappedReads) / transcript.RefLength;
        }

        // Now posterior has the transcript fraction
        for (auto& transcript : transcripts_) {
            double logLength = std::log(transcript.RefLength);
            double fpkmFactor = std::exp(logBillion - logLength - logNumFragments);
            double count = transcript.projectedCounts;
            //double countTotal = transcripts_[transcriptID].totalCounts;
            //double countUnique = transcripts_[transcriptID].uniqueCounts;
            double fpkm = count > 0 ? fpkmFactor * count : 0.0;
            double npm = (transcript.projectedCounts / numMappedReads);
            double tfrac = (npm / transcript.RefLength) / tfracDenom;
            double tpm = tfrac * million;

            fmt::print(output.get(), "{}\t{}\t{}\t{}\t{}\n",
                    transcript.RefName, transcript.RefLength,
                    tpm, fpkm, count);
        }

    }

    LibraryFormat hitType(uint32_t end1Start, bool end1Fwd,
                          uint32_t end2Start, bool end2Fwd) {

        // If the reads come from opposite strands
        if (end1Fwd != end2Fwd) {
            // and if read 1 comes from the forward strand
            if (end1Fwd) {
                // then if read 1 start < read 2 start ==> ISF
                if (end1Start <= end2Start) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::SA);
                } // otherwise read 2 start < read 1 start ==> OSF
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::SA);
                }
            }
            // and if read 2 comes from the forward strand
            if (end2Fwd) {
                // then if read 2 start <= read 1 start ==> ISR
                if (end2Start <= end1Start) {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD, ReadStrandedness::AS);
                } // otherwise, read 2 start > read 1 start ==> OSR
                else {
                    return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::AWAY, ReadStrandedness::AS);
                }
            }
        } else { // Otherwise, the reads come from the same strand
            if (end1Fwd) { // if it's the forward strand ==> MSF
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::S);
            } else { // if it's the reverse strand ==> MSR
                return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::SAME, ReadStrandedness::A);
            }
        }
        // SHOULD NOT GET HERE
        fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);
    }


    LibraryFormat hitType(uint32_t start, bool isForward) {
        // If the read comes from the forward strand
        if (isForward) {
            return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S);
        } else {
            return LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE, ReadStrandedness::S);
        }
        // SHOULD NOT GET HERE
        fmt::print(stderr, "WARNING: Could not associate known library type with read!\n");
        return LibraryFormat(ReadType::PAIRED_END, ReadOrientation::NONE, ReadStrandedness::U);

    }

}
}

template
void salmon::utils::writeAbundances<AlignmentLibrary<ReadPair>>(AlignmentLibrary<ReadPair>& alnLib,
                                              boost::filesystem::path& fname,
                                              std::string headerComments);

template
void salmon::utils::writeAbundances<AlignmentLibrary<UnpairedRead>>(AlignmentLibrary<UnpairedRead>& alnLib,
                                                  boost::filesystem::path& fname,
                                                  std::string headerComments);
template
void salmon::utils::writeAbundances<ReadExperiment>(ReadExperiment& alnLib,
                                                  boost::filesystem::path& fname,
                                                  std::string headerComments);

