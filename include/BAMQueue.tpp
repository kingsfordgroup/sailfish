#include "BAMQueue.hpp"

template <typename FragT>
BAMQueue<FragT>::BAMQueue(std::vector<boost::filesystem::path>& fnames, LibraryFormat& libFmt):
    files_(std::vector<AlignmentFile>()),
    libFmt_(libFmt), totalReads_(0),
    numUnaligned_(0), numMappedReads_(0), doneParsing_(false) {

        size_t capacity{2000000};
        fragmentQueue_.set_capacity(capacity);
        for (size_t i = 0; i < capacity; ++i) {
            fragmentQueue_.push(new FragT);
        }

        size_t groupCapacity = 5000000;
        alnGroupPool_.set_capacity(groupCapacity);
        for (size_t i = 0; i < groupCapacity; ++i) {
            alnGroupPool_.push(new AlignmentGroup<FragT*>);
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
  fmt::print(stderr, "Resetting BAMQueue from file [{}] . . .", fname_);
  parsingThread_->join();
  for (auto& file : files_) {
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
  totalReads_ = 0;
  numUnaligned_ = 0;
  numMappedReads_ = 0;
  doneParsing_ = false;
  batchNum_ = 0;
  fmt::print(stderr, "done\n");
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
    while (!doneParsing_ or !alnGroupQueue_.empty()) {
        if (alnGroupQueue_.pop(group)) {
            return true;
        }
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
void BAMQueue<FragT>::start() {
    fmt::print(stderr, "Started parsing\n");
    // Depending on the specified library type, start an
    // appropriate parsing thread.
    parsingThread_.reset(new std::thread([this]()-> void {
            this->fillQueue_();
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

template <typename FragT>
inline bool BAMQueue<FragT>::getFrag_(ReadPair& rpair) {
    bool haveValidPair{false};

       while (!haveValidPair) {
        bool didRead1 = (sam_read1(fp_, hdr_, rpair.read1) >= 0);
        bool didRead2 = (sam_read1(fp_, hdr_, rpair.read2) >= 0);
        // If we didn't get a read, then we've exhausted this file
        // NOTE: Not sure about the *or* condition here. In some cases,
        // we may be discarding a single read, but it won't be properly
        // paired anyway.  Figure out what the right thing is to do 
        // here.
        if (!didRead1 or !didRead2) { 
            currFile_++;
            // If this is the last file, then we're done
            if (currFile_ == files_.end()) { return false; }
            // Otherwise, start parsing the next file.
            fp_ = currFile_->fp;
            hdr_ = currFile_->header;
            continue;
        }

        bool read1IsValid{false};
        if (rpair.read1->core.flag & BAM_FPROPER_PAIR 
            and !(rpair.read1->core.flag & BAM_FDUP) 
            and !(rpair.read1->core.flag & BAM_FQCFAIL)
            ) {
            read1IsValid = true;
        }

        bool read2IsValid{false};
        if (rpair.read2->core.flag & BAM_FPROPER_PAIR 
            and !(rpair.read2->core.flag & BAM_FDUP) 
            and !(rpair.read2->core.flag & BAM_FQCFAIL)
            ) {
            read2IsValid = true;
        }

        haveValidPair = read1IsValid and read2IsValid and
            (rpair.read1->core.tid == rpair.read2->core.tid) and
            (rpair.read1->core.l_qname == rpair.read2->core.l_qname) and
            (strcmp(bam1_qname(rpair.read1), bam1_qname(rpair.read2)) == 0);

        if (!haveValidPair) { ++numUnaligned_; }//bam_destroy1(rpair.read1); bam_destroy1(rpair.read2); ++numUnaligned_;}
        ++totalReads_;
    }

    rpair.logProb = sailfish::math::LOG_0;
    return true;
}

template <typename FragT>
inline bool BAMQueue<FragT>::getFrag_(UnpairedRead& sread) {
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
            sread.transcriptID() >= 0) {
            haveValidRead = true;
        }

        if (!haveValidRead) { ++numUnaligned_; }
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
void BAMQueue<FragT>::fillQueue_() {
    size_t n{0};
    //AlignmentGroup* alngroup = new AlignmentGroup;
    AlignmentGroup<FragT*>* alngroup;
    alnGroupPool_.pop(alngroup);
     
    currFile_ = files_.begin();
    fp_ = currFile_->fp;
    hdr_ = currFile_->header;

    FragT* f;
    fragmentQueue_.pop(f);
    //ReadPair p = {nullptr, nullptr, LOG_0};
    char* prevReadName = new char[100];
    prevReadName[0] = '\0';
    while(getFrag_(*f) and !doneParsing_) {

        char* readName = f->getName();//bam1_qname(p.read1);
        // if this is a new read
        // TODO: (p.read1->core.l_qname != p.read2->core.l_qname)
        if (strcmp(readName, prevReadName) != 0) {
            if (alngroup->size() > 0)  {
                // push the align group
                while(!alnGroupQueue_.push(alngroup));
                alngroup = nullptr;
                alnGroupPool_.pop(alngroup);
            }
            alngroup->addAlignment(f);
            prevReadName = readName;
            numMappedReads_++;
        } else {
            alngroup->addAlignment(f);
       }
        fragmentQueue_.pop(f);
        ++n;
    }
    currFile_ = files_.end();
    fp_ = nullptr;
    hdr_ = nullptr;
    doneParsing_ = true;
    return;
}


