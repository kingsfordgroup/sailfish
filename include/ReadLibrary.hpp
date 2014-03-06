#ifndef READ_LIBRARY_HPP
#define READ_LIBRARY_HPP

#include <vector>
#include <exception>

#include "LibraryFormat.hpp"

/**
 * This class represents the basic information about a library of reads, like
 * its paired-end status, the reads that should appear on the forward and reverse strand,
 * and the relative orientation of the reads. 
 */
class ReadLibrary {
public:
    /**
     * Construct a new ReadLibrary of the given format
     */
    ReadLibrary(LibraryFormat& fmt) : fmt_(fmt) {}

    /**
     * Add files containing mated reads (from pair 1 of the mates) to this library. 
     */
    void addMates1(const std::vector<std::string>& mateOneFilenames) {
        mateOneFilenames_ = mateOneFilenames;
    }
    
    /**
     * Add files containing mated reads (from pair 2 of the mates) to this library.
     */
    void addMates2(const std::vector<std::string>& mateTwoFilenames) {
        mateTwoFilenames_ = mateTwoFilenames;
    }
    
    /**
     * Add files containing unmated reads.
     */
    void addUnmated(const std::vector<std::string>& unmatedFilenames) {
        unmatedFilenames_ = unmatedFilenames;
    }

    /**
     * Return true if this read library is for paired-end reads and false otherwise.
     */
    bool isPairedEnd() {
        return (fmt_.type == ReadType::PAIRED_END); 
    }

    /**
     * Checks if this read library is valid --- if it's paired-end, it should have mate1/2 reads and the same
     * number of files for each; if it's unpaired it should have only unpaired files.  
     * NOTE: This function throws an exception if this is not a valid read library!
     */
    void checkValid() {
        if (isPairedEnd()) {
            size_t n1 = mateOneFilenames_.size();
            size_t n2 = mateTwoFilenames_.size();
            if (n1 == 0 or n2 == 0 or n1 != n2) {
                std::string e = "You must provide #1 and #2 mated read files with a paired-end library type";
                throw std::invalid_argument(e);
            }
        } else {
            size_t n = unmatedFilenames_.size();
            if (n == 0) {
                std::string e= "You must provide unmated read files with a single-end library type";
                throw std::invalid_argument(e);
            }

        }
    }

    /**
     * Return the vector of files containing the mate1 reads for this library.
     */
    const std::vector<std::string>& mates1() const { return mateOneFilenames_; }

    /**
     * Return the vector of files containing the mate2 reads for this library.
     */
    const std::vector<std::string>& mates2() const { return mateTwoFilenames_; }

    /**
     * Return the vector of files containing the unmated reads for the library.
     */
    const std::vector<std::string>& unmated() const { return unmatedFilenames_; }
    
    /**
     * Return the LibraryFormat object describing the format of this read library.
     */
    const LibraryFormat& format() const { return fmt_; }
    
private:
    LibraryFormat fmt_;
    std::vector<std::string> unmatedFilenames_;
    std::vector<std::string> mateOneFilenames_;
    std::vector<std::string> mateTwoFilenames_;
    bool isPairedEnd_;
};

#endif // READ_LIBRARY_HPP
