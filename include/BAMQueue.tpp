#include "BAMQueue.hpp"

template <typename FragT>
BAMQueue<FragT>::BAMQueue(const std::string& fname, LibraryFormat& libFmt):
    fname_(fname), libFmt_(libFmt), totalReads_(0),
    numUnaligned_(0), doneParsing_(false) {

        size_t capacity{10000000};
        alnStructQueue_.set_capacity(capacity);
        for (size_t i = 0; i < capacity; ++i) {
            alnStructQueue_.push(bam_init1());
        }

        size_t groupCapacity = 5000000;
        alnGroupPool_.set_capacity(groupCapacity);
        for (size_t i = 0; i < groupCapacity; ++i) {
            alnGroupPool_.push(new AlignmentGroup<FragT>);
        }


        auto hasSuffix = [](const std::string &str, const std::string &suffix) -> bool {
                return str.size() >= suffix.size() &&
                str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
            };

        std::string readMode = "r";
        if (hasSuffix(fname, ".bam")) {
            readMode = "rb";
        }
        //bam_set_num_threads_per(8);
        fp_ = hts_open(fname_.c_str(), readMode.c_str());
        hts_set_threads(fp_, 8);
        hdr_ = sam_hdr_read(fp_);
}

template <typename FragT>
BAMQueue<FragT>::~BAMQueue() {
    parsingThread_->join();
    sam_close(fp_);
    //samclose(fp_);

    // Free the structure holding all of the reads
    bam1_t* aln;
    while(!alnStructQueue_.empty()) { alnStructQueue_.pop(aln); bam_destroy1(aln); }
    AlignmentGroup<FragT>* grp;
    while(!alnGroupPool_.empty()) { alnGroupPool_.pop(grp); delete grp; }

}

template <typename FragT>
inline bool BAMQueue<FragT>::getAlignmentGroup(AlignmentGroup<FragT>*& group) {
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
bam_header_t* BAMQueue<FragT>::header() { return hdr_; } //fp_->header; }

template <typename FragT>
void BAMQueue<FragT>::start() {
    // Depending on the specified library type, start an
    // appropriate parsing thread.
    parsingThread_ = new std::thread([this]()-> void {
            this->fillQueue_();
    });
}

template <typename FragT>
tbb::concurrent_bounded_queue<bam1_t*>& BAMQueue<FragT>::getAlignmentStructureQueue() {
    return alnStructQueue_;
}

template <typename FragT>
tbb::concurrent_bounded_queue<AlignmentGroup<FragT>*>& BAMQueue<FragT>::getAlignmentGroupQueue() {
    return alnGroupPool_;
}

template <typename FragT>
inline bool BAMQueue<FragT>::getFrag_(ReadPair& rpair) {
    bool haveValidPair{false};
    alnStructQueue_.pop(rpair.read1);
    alnStructQueue_.pop(rpair.read2);

    while (!haveValidPair) {

        bool didRead1 = (sam_read1(fp_, hdr_, rpair.read1) >= 0);
        if (!didRead1) { alnStructQueue_.push(rpair.read1); alnStructQueue_.push(rpair.read2); return false; }

        bool didRead2 = (sam_read1(fp_, hdr_, rpair.read2) >= 0);
        if (!didRead2) { alnStructQueue_.push(rpair.read1); alnStructQueue_.push(rpair.read2); return false; }

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
    alnStructQueue_.pop(sread.read);

    while (!haveValidRead) {
        bool didRead = (sam_read1(fp_, hdr_, sread.read) >= 0);//(samread(fp_, sread.read) > 0);
        if (!didRead) { alnStructQueue_.push(sread.read); return false; }

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
void BAMQueue<FragT>::fillQueue_() {
    size_t n{0};
    //AlignmentGroup* alngroup = new AlignmentGroup;
    AlignmentGroup<FragT>* alngroup;
    alnGroupPool_.pop(alngroup);
    FragT f;
    //ReadPair p = {nullptr, nullptr, LOG_0};
    char* prevReadName = new char[100];
    prevReadName[0] = '\0';
    while(getFrag_(f) and !doneParsing_) {

        char* readName = f.getName();//bam1_qname(p.read1);
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
        } else {
            alngroup->addAlignment(f);
        }
        ++n;
    }
    doneParsing_ = true;
    return;
}


