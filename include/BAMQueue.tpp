#include "BAMQueue.hpp"
#include "IOUtils.hpp"
#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY
#include <chrono>

#include "g2logworker.h"
#include "g2log.h"

template <typename FragT>
BAMQueue<FragT>::BAMQueue(std::vector<boost::filesystem::path>& fnames, LibraryFormat& libFmt):
    files_(std::vector<AlignmentFile>()),
    libFmt_(libFmt), totalReads_(0),
    numUnaligned_(0), numMappedReads_(0), doneParsing_(false) {

        size_t capacity{2000000};
        fragmentQueue_.set_capacity(capacity);
        for (size_t i = 0; i < capacity; ++i) {
            // avoid r-value ref until we figure out what's
            // up with TBB 4.3
            auto* fragPtr = new FragT;
            fragmentQueue_.push(fragPtr);
        }

        size_t groupCapacity = 5000000;
        alnGroupPool_.set_capacity(groupCapacity);
        for (size_t i = 0; i < groupCapacity; ++i) {
            // avoid r-value ref until we figure out what's
            // up with TBB 4.3
            auto* agrpPtr = new AlignmentGroup<FragT*>; 
            alnGroupPool_.push(agrpPtr);
        }

        for (auto& fname : fnames) {
            readMode_ = "r";
            if (fname.extension() == ".bam") {
                readMode_ = "rb";
            }
            auto* fp = hts_open(fname.c_str(), readMode_.c_str());
            files_.push_back({fname, fp, sam_hdr_read(fp)});
        }
}

template <typename FragT>
void BAMQueue<FragT>::reset() {
  fmt::print(stderr, "Resetting BAMQueue from file(s) [ ");
  parsingThread_->join();
  for (auto& file : files_) {
    fmt::print(stderr, "{} ", file.fileName);
    // re-open the file
    hts_close(file.fp); 
    file.fp = nullptr;
    readMode_ = (file.fileName.extension() == ".bam") ? "rb" : "r";
    file.fp = hts_open(file.fileName.c_str(), readMode_.c_str());
    // re-obtain the header
    bam_header_destroy(file.header); 
    file.header = nullptr;
    file.header = sam_hdr_read(file.fp);
  }
  fmt::print(stderr, "] . . . done\n");
  totalReads_ = 0;
  numUnaligned_ = 0;
  numMappedReads_ = 0;
  doneParsing_ = false;
  batchNum_ = 0;
}

template <typename FragT>
BAMQueue<FragT>::~BAMQueue() {
    fmt::print(stderr, "\nFreeing memory used by read queue . . . ");
    parsingThread_->join();
    for (auto& file : files_) {
        // close the file 
        hts_close(file.fp); 
        file.fp = nullptr;
        // destroy the header
        bam_header_destroy(file.header); 
        file.header= nullptr;
    }

    // Free the structure holding all of the reads
    FragT* frag;
    while(!fragmentQueue_.empty()) { 
        fragmentQueue_.pop(frag); 
        delete frag;
        frag = nullptr;
    }

    AlignmentGroup<FragT*>* grp;
    while(!alnGroupPool_.empty()) { alnGroupPool_.pop(grp); delete grp; grp = nullptr; }
    while(!alnGroupQueue_.empty()) { alnGroupQueue_.pop(grp); delete grp; grp = nullptr; }
    fmt::print(stderr, "done\n");
}

template <typename FragT>
inline bool BAMQueue<FragT>::getAlignmentGroup(AlignmentGroup<FragT*>*& group) {
    volatile bool isEmpty{alnGroupQueue_.empty()};
    while (!doneParsing_) {
        while (alnGroupQueue_.pop(group)) {
            return true;
        }
    }
    while (alnGroupQueue_.pop(group)) {
        return true;
    }
    return false;
}


template <typename FragT>
void BAMQueue<FragT>::forceEndParsing() { doneParsing_ = true; }

template <typename FragT>
bam_header_t* BAMQueue<FragT>::header() { return files_.front().header; } 

template <typename FragT>
std::vector<bam_header_t*> BAMQueue<FragT>::headers() { 
    std::vector<bam_header_t*> hs;
    for (auto& file : files_) {
        hs.push_back(file.header);
    }
    return hs;
}

template <typename FragT>
template <typename FilterT>
void BAMQueue<FragT>::start(FilterT filt) {
    fmt::print(stderr, "Started parsing\n");
    // Depending on the specified library type, start an
    // appropriate parsing thread.
    parsingThread_.reset(new std::thread([this, filt]()-> void {
            this->fillQueue_(filt);
    }));
}

template <typename FragT>
tbb::concurrent_bounded_queue<FragT*>& BAMQueue<FragT>::getFragmentQueue() {
    return fragmentQueue_;
}

template <typename FragT>
tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>& BAMQueue<FragT>::getAlignmentGroupQueue() {
    return alnGroupPool_;
}

inline bool checkProperPairedNames_(const char* qname1, const char* qname2, const uint32_t nameLen) {
    bool same{true};
    bool sameEnd{true};
    if (BOOST_LIKELY(nameLen > 1)) {
        same = (memcmp(qname1, qname2, nameLen - 1) == 0);
        sameEnd = qname1[nameLen-1] == qname2[nameLen-1];
        same = same and (sameEnd or qname1[nameLen-1] == '1' or qname1[nameLen-1] == '2');
        same = same and (sameEnd or qname2[nameLen-1] == '2' or qname2[nameLen-1] == '1');
    } else {
        same = (memcmp(qname1, qname2, nameLen) == 0);
    }
    return same;
}


/**
* This set of functions checks if two reads have the same name.  If the reads
* are paired-end, it allows them to differ in the last character if they end in 1/2.
* If the reads are single-end, the strings must be identical.
*/
template <typename T>
inline bool sameReadName_(const char* qname1, const char* qname2, const uint32_t nameLen)  = delete;


// Specialization for unpaired reads.
template <>
inline bool sameReadName_<UnpairedRead>(const char* qname1, const char* qname2, const uint32_t nameLen) {
    return (memcmp(qname1, qname2, nameLen) == 0);
}

// Specialization for paired reads.
template <>
inline bool sameReadName_<ReadPair>(const char* qname1, const char* qname2, const uint32_t nameLen) {
    return checkProperPairedNames_(qname1, qname2, nameLen);
}

enum class AlignmentType : uint8_t { 
    UnmappedOrphan = 0,
    MappedOrphan = 1,
    MappedDiscordantPair = 2,
    MappedConcordantPair = 3,
    UnmappedPair = 4
};

inline AlignmentType getPairedAlignmentType_(bam1_t* aln) {
    bool readIsMapped = !(aln->core.flag & BAM_FUNMAP);
    bool mateIsMapped = !(aln->core.flag & BAM_FMUNMAP);

    if (readIsMapped and mateIsMapped) {
        if (aln->core.flag & BAM_FPROPER_PAIR) {
            return AlignmentType::MappedConcordantPair;
        } else {
            return AlignmentType::MappedDiscordantPair;
        }
    } 
    if (readIsMapped and !mateIsMapped) {
        return AlignmentType::MappedOrphan;
    }
    if (mateIsMapped and !readIsMapped) {
        return AlignmentType::UnmappedOrphan;
    }
    if (!mateIsMapped and !readIsMapped) {
        return AlignmentType::UnmappedPair;
    }
    std::cerr << "\n\n\nEncountered unknown alignemnt type; this should not happen!\n"
              << "Please file a bug report on GitHub. Exiting.\n";
    std::exit(1);
}

template <typename FragT>
template <typename FilterT>
inline bool BAMQueue<FragT>::getFrag_(ReadPair& rpair, FilterT filt) {
    bool haveValidPair{false};
    bool didRead1{false};
    bool didRead2{false};

    rpair.orphanStatus = salmon::utils::OrphanStatus::LeftOrphan;
    while (!haveValidPair) {
        // Consume a single read
        didRead1 = (sam_read1(fp_, hdr_, rpair.read1) >= 0);

        AlignmentType alnType;
        // If we were able to obtain a read, determine what type 
        // of alignment it came from.
        while (didRead1) {
            alnType = getPairedAlignmentType_(rpair.read1);
            // If this read is part of a concordantly mapped pair, or an
            // unmapped pair, then go get the other read.
            if (   (alnType == AlignmentType::MappedConcordantPair)
                or (alnType == AlignmentType::UnmappedPair)) { break; }

            switch (alnType) {
                case AlignmentType::UnmappedOrphan:
                    ++numUnaligned_;
                    if (filt != nullptr) {
                        rpair.orphanStatus = salmon::utils::OrphanStatus::LeftOrphan;
                        filt->processFrag(&rpair);
                    }
                    break;
                // === end of UnmappedOrphan case
                case AlignmentType::MappedOrphan:
                    rpair.orphanStatus =  (rpair.read1->core.flag & BAM_FREVERSE) ? 
                                          salmon::utils::OrphanStatus::LeftOrphan :
                                          salmon::utils::OrphanStatus::RightOrphan;
                    rpair.logProb = sailfish::math::LOG_0;
                    return true;
                // === end of MappedOrphan case 
                case AlignmentType::MappedDiscordantPair:
                    // The discordant mapped pair case is sort of a nightmare
                    // it requires arbitrary look-ahead followed by a pairing 
                    // of the ends once all of the records for a read have been
                    // consumed.  For the time being, we don't support discordant
                    // mappings --- just skip over them.

                    // Don't even pass them to the filter right now!
                    /*
                    if (filt != nullptr) {
                        rpair.orphanStatus = salmon::utils::OrphanStatus::LeftOrphan;
                        filt->processFrag(&rpair);
                    }
                    */
                    break;
                // === end if MappedDiscordantPair case
                default:
                    std::cerr << "\n\n\nEncountered unknown alignemnt type; this should not happen!\n"
                        << "Please file a bug report on GitHub. Exiting.\n";
                    std::exit(1);
                    break;
            }
            // If this was not a properly mapped orphan read, then grab the next
            // read.
            didRead1 = (sam_read1(fp_, hdr_, rpair.read1) >= 0);
        }

        didRead2 = (sam_read1(fp_, hdr_, rpair.read2) >= 0);

        // If we didn't get a read, then we've exhausted this file. 
        // NOTE: I'm not sure about the *or* condition here. In some cases, we
        // may be discarding a single read, but it won't be properly paired
        // anyway. Figure out what the right thing is to do here.
        if (!didRead1 or !didRead2) { 
            currFile_++;
            // If this is the last file, then we're done
            if (currFile_ == files_.end()) { return false; }
            // Otherwise, start parsing the next file.
            fp_ = currFile_->fp;
            hdr_ = currFile_->header;
            continue;
        }

        
        // If we expected a paired read, but didn't find one
        // then flip out and quit!
        if (BOOST_UNLIKELY((
            !(rpair.read1->core.flag & BAM_FPAIRED) and
            !(rpair.read2->core.flag & BAM_FPAIRED))
            )) {
            fmt::MemoryWriter errmsg;
            errmsg << "\n\n" 
                   << ioutils::SET_RED << "ERROR: " << ioutils::RESET_COLOR 
                   << "Found unpaired read in a paired-end library. "
                   << "The read was marked as unpaired in sequencing (not just unmapped)." 
                   << "The two ends of a paired-end read should be adjacent. "
                   << "Don't know how to proceed; exiting!\n\n";
            fmt::print(stderr, errmsg.str());
            LOG(WARNING) << errmsg.str();
            std::exit(-1);
        }
        
        // We've observed two, consecutive paired reads; now check if our reads
        // have the same name.

        // The names must first have the same length.
        bool sameName = (rpair.read1->core.l_qname == rpair.read2->core.l_qname);

        // If the lengths are the same, check the actual strings.  Use
        // memcmp for efficiency since we know the length.
        if (BOOST_LIKELY(sameName)) {
            auto nameLen = rpair.read1->core.l_qname;
            char* qname1 = bam1_qname(rpair.read1);
            char* qname2 = bam1_qname(rpair.read2);
            sameName = (memcmp(qname1, qname2, nameLen) == 0);
        }

        // If we've gotten this far, then the read should be a pair (same name)
        // and should either be concordantly aligned (proper pair) or both unaligned.
        if (BOOST_UNLIKELY(
                !sameName or ((rpair.read1->core.flag & BAM_FPROPER_PAIR) != (rpair.read2->core.flag & BAM_FPROPER_PAIR)))) {
            std::cerr << "\n\n\n";
            std::cerr << "WARNING: Detected suspicious pair --- \n";
            if (!sameName) {
                std::cerr << "\tThe names are different:\n";
                std::cerr << "\tread1 : " << bam1_qname(rpair.read1) << "\n";
                std::cerr << "\tread2 : " << bam1_qname(rpair.read2) << "\n";
            }
            if((rpair.read1->core.flag & BAM_FPROPER_PAIR) != (rpair.read2->core.flag & BAM_FPROPER_PAIR)) {
                std::cerr << "\tThe proper-pair statuses are inconsistent:\n"; 
                std::cerr << "read1 : " << (!(rpair.read1->core.flag & BAM_FPROPER_PAIR) ? "no " : "") << "proper-pair; "
                          << ((rpair.read1->core.flag & BAM_FUNMAP) ? "not " : "") << "mapped; mate"
                          << ((rpair.read1->core.flag & BAM_FMUNMAP) ? "not " : "") << "mapped\n\n";
                std::cerr << "read2: " << (!(rpair.read2->core.flag & BAM_FPROPER_PAIR) ? "no " : "") << "PROPER PAIR; "
                          << ((rpair.read2->core.flag & BAM_FUNMAP) ? "not " : "") << "mapped; mate"
                          << ((rpair.read2->core.flag & BAM_FMUNMAP) ? "not " : "") << "mapped\n\n";
            }
        }

        switch (alnType) {
            case AlignmentType::MappedConcordantPair:
                haveValidPair = true;
                if (rpair.read1->core.flag & BAM_FREAD2) {
                    std::swap(rpair.read1, rpair.read2);
                }
                rpair.orphanStatus = salmon::utils::OrphanStatus::Paired;
                break;
            case AlignmentType::UnmappedPair:
                ++numUnaligned_;
                if ((filt != nullptr) and sameName) {
                    rpair.orphanStatus = salmon::utils::OrphanStatus::Paired;
                    filt->processFrag(&rpair);
                }
                break;
            default:
                std::cerr << "\n\n\nEncountered unknown alignemnt type; this should not happen!\n"
                          << "Please file a bug report on GitHub. Exiting.\n";
                std::exit(1);
                break;
        }
        ++totalReads_;
    }
    rpair.logProb = sailfish::math::LOG_0;
    return true;
}

template <typename FragT>
template <typename FilterT>
inline bool BAMQueue<FragT>::getFrag_(UnpairedRead& sread, FilterT filt) {
    bool haveValidRead{false};

    while (!haveValidRead) {
        bool didRead = (sam_read1(fp_, hdr_, sread.read) >= 0);
        // If we didn't get a read, then we've exhausted this file
        if (!didRead) { 
            currFile_++;
            // If this is the last file, then we're done
            if (currFile_ == files_.end()) { return false; }
            // Otherwise, start parsing the next file.
            fp_ = currFile_->fp;
            hdr_ = currFile_->header;
        }

        if (!(sread.read->core.flag & BAM_FDUP) and
            !(sread.read->core.flag & BAM_FQCFAIL) and
            !(sread.read->core.flag & BAM_FUNMAP) and
            sread.transcriptID() >= 0) {
            haveValidRead = true;
        }

        if (!haveValidRead) { 
            if (filt != nullptr) {
                filt->processFrag(&sread);
            }
            ++numUnaligned_; 
        }
        ++totalReads_;
    }

    sread.logProb = sailfish::math::LOG_0;
    return true;
}

template <typename FragT>
size_t BAMQueue<FragT>::numObservedReads(){ return totalReads_; }

template <typename FragT>
size_t BAMQueue<FragT>::numMappedReads(){ 
    return numMappedReads_;
}

template <typename FragT>
template <typename FilterT>
void BAMQueue<FragT>::fillQueue_(FilterT filt) {
    size_t n{0};
    AlignmentGroup<FragT*>* alngroup;
    alnGroupPool_.pop(alngroup);

    currFile_ = files_.begin();
    fp_ = currFile_->fp;
    hdr_ = currFile_->header;

    FragT* f;
    fragmentQueue_.pop(f);

    uint32_t prevLen{1};
    char* prevReadName = new char[255];
    
    while(getFrag_(*f, filt)) {
        char* readName = f->getName();
        uint32_t currLen = f->getNameLength();
        std::string readStr(readName);
        // if this is a new read
        if ( (currLen != prevLen) or
            !sameReadName_<UnpairedRead>(readName, prevReadName, currLen) ) {

            if (alngroup->size() > 0) {
                // push the align group
                while(!alnGroupQueue_.push(alngroup));
                alngroup = nullptr;
                alnGroupPool_.pop(alngroup);
            }
            alngroup->addAlignment(f);
            f = nullptr;
            memcpy(prevReadName, readName, currLen);
            prevLen = currLen;
            numMappedReads_++;
        } else { // otherwise, this is another alignment for the same read
            alngroup->addAlignment(f);
            f = nullptr;
        }
        fragmentQueue_.pop(f);
        ++n;
    }

    // If we popped a fragment structure off the queue, but didn't add it
    // to an alignment group, then reclaim it here
    if (f != nullptr) { fragmentQueue_.push(f); f = nullptr; }

    // If the last alignment group is non-empty, then send
    // it off to be processed.
    if (alngroup->size() > 0) {
        while(!alnGroupQueue_.push(alngroup));
        alngroup = nullptr;
    } else { // otherwise, reclaim the alignment group structure here
        alnGroupPool_.push(alngroup);
    }

    delete [] prevReadName;

    // We're at the end of the list of input files
    // and we're done parsing (for now).
    currFile_ = files_.end();
    fp_ = nullptr;
    hdr_ = nullptr;
    doneParsing_ = true;
    return;
}

///////// Proper BAM parsing graveyard

/* 
        // If the reads don't have the same name, then the pair was not
        // consecutive in the file --- complain and skip!
        if (BOOST_UNLIKELY(!sameName)) {

            // Flag stores whether or not we've consumed all of the 
            // alignment records in the current file.
            bool consumedFile{false};

            // Consume reads until we find a pair with the same name
            while (!sameName) {
                // Complain if this is supposed to be a paired read
                fmt::print(stderr, "{}WARNING:{} The mate of read [{}] did not "
                        "appear next to it in the file. The next read was [{}].  "
                        "Skipping the first read\n\n", 
                        ioutils::SET_RED, 
                        ioutils::RESET_COLOR, 
                        bam1_qname(rpair.read1),
                        bam1_qname(rpair.read2));

                rpair.read1 = rpair.read2; 
                didRead2 = (sam_read1(fp_, hdr_, rpair.read2) >= 0);
                // If we hit the end of the file --- skip to the top of the loop
                // to see if we need to move on to another file.
                if (!didRead2) { consumedFile = true; continue; }

                // As above, if we encounter a non-paired read, then, for the
                // time being, flip out and quit.
                if (BOOST_UNLIKELY(!(rpair.read2->core.flag & BAM_FPAIRED))) {
                    fmt::MemoryWriter errmsg;
                    errmsg << "\n\n" 
                        << ioutils::SET_RED << "ERROR: " << ioutils::RESET_COLOR 
                        << "Saw adjacent reads, at least one of which was unpaired. "
                        << "The two ends of a paired-end read should be adjacent. "
                        << "Don't know how to proceed; exiting!\n\n";
                    std::cerr << errmsg.str();
                    std::exit(-1);
                } else { std::cerr << "BLARGHHH\n\n\n"; }
                // As above, check first that the lengths of the names are the
                // same and then that the names are, in fact, identical.
                sameName = (rpair.read1->core.l_qname == rpair.read2->core.l_qname);
                if (BOOST_LIKELY(sameName)) {
                    auto nameLen = rpair.read1->core.l_qname;
                    char* qname1 = bam1_qname(rpair.read1);
                    char* qname2 = bam1_qname(rpair.read2);
                    sameName = checkProperPairedNames_(qname1, qname2, nameLen);
                }
            } // end while (!sameName)

            // If we consumed all of the current file, break to the top of the
            // loop.
            if (BOOST_UNLIKELY(consumedFile)) { continue; } 
        }

        bool read1IsValid{false};
        if ( !(rpair.read1->core.flag & BAM_FUNMAP) 
              and (rpair.read1->core.flag & BAM_FPROPER_PAIR)
             //and !(rpair.read1->core.flag & BAM_FDUP) 
             //and !(rpair.read1->core.flag & BAM_FQCFAIL)
            ) {
            read1IsValid = true;
        }

        bool read2IsValid{false};
        if ( !(rpair.read2->core.flag & BAM_FUNMAP) 
              and (rpair.read2->core.flag & BAM_FPROPER_PAIR)
             //and !(rpair.read2->core.flag & BAM_FDUP) 
             //and !(rpair.read2->core.flag & BAM_FQCFAIL)
            ) {
            read2IsValid = true;
        }

        haveValidPair = read1IsValid and read2IsValid and
                        (rpair.read1->core.tid == rpair.read2->core.tid) and 
                        sameName;
        
        // If the pair was not properly mapped 
        if (!haveValidPair) { 
            ++numUnaligned_; 
            // If we have an active output filter, and the reads were not 
            // a valid alignment (*but were a proper pair*), then pass them
            // to the output filter.
            if ((filt != nullptr) and sameName) {
                filt->processFrag(&rpair);
            }
        } else {
            // Make sure read1 is read1 and read2 is read2; else swap
            if (rpair.read1->core.flag & BAM_FREAD2) {
                std::swap(rpair.read1, rpair.read2);
            }
            rpair.orphanStatus = salmon::utils::OrphanStatus::Paired;

            ++propPaired;
        }
        */



