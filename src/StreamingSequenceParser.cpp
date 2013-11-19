/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/


#include "StreamingSequenceParser.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>
#include "unistd.h"
#include "fcntl.h"
#include <poll.h>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include "tbb/concurrent_queue.h"


int fifoRead(int f, void* buf, int n) {
  char* a;
  int m, t;
  a = static_cast<char*>(buf);
  t = 0;
  while (t < n) {
    m = read(f, a+t, n-t);
    if (m <= 0) {
      if (t == 0) {
        return m;
      }
      std::cerr << "hit eof!\n";
      break;
    }
    t += m;
  }
  return t;
}

KSEQ_INIT(int, fifoRead)

namespace bfs = boost::filesystem;

StreamingReadParser::StreamingReadParser( std::vector<bfs::path>& files ): inputStreams_(files),
        parsing_(false), parsingThread_(nullptr)
    {
        readStructs_ = new ReadSeq[queueCapacity_];
        readQueue_.set_capacity(queueCapacity_);
        seqContainerQueue_.set_capacity(queueCapacity_);
        for (size_t i = 0; i < queueCapacity_; ++i) {
            seqContainerQueue_.push(&readStructs_[i]);
        }
    }

StreamingReadParser::~StreamingReadParser() {
        parsingThread_->join();
        for (auto i : boost::irange(size_t{0}, queueCapacity_)) {
            if (readStructs_[i].seq != nullptr) { free(readStructs_[i].seq); }
        }
        delete [] readStructs_;
        delete parsingThread_;
    }

bool StreamingReadParser::start() {
        if (!parsing_) {
            parsing_ = true;
            parsingThread_ = new std::thread([this](){

                kseq_t* seq;
                ReadSeq* s;
                std::cerr << "reading from " << this->inputStreams_.size() << " streams\n";
                for (auto file : this->inputStreams_) {
                    std::cerr << "reading from " << file.native() << "\n";
                    // open the file and init the parser
                    struct pollfd pfd;
                    pfd.events = POLLERR;
                    pfd.fd = open(file.c_str(), O_RDONLY);


                    //FILE* fdstream = fdopen(fp, "r");
                    seq = kseq_init(pfd.fd);
                      int ksv = kseq_read(seq);
                      while (ksv >= 0) {
                        this->seqContainerQueue_.pop(s);

                        // Possibly allocate more space for the sequence
                        if (seq->seq.l > s->len) {
                            s->seq = static_cast<char*>(realloc(s->seq, seq->seq.l));
                        }
                        // Copy the sequence length and sequence over to the ReadSeq struct
                        s->len = seq->seq.l;
                        memcpy(s->seq, seq->seq.s, s->len);

                        // Possibly allocate more space for the name
                        if (seq->name.l > s->nlen) {
                            s->name = static_cast<char*>(realloc(s->name, seq->name.l));
                        }
                        // Copy the name length and name over to the ReadSeq struct
                        s->nlen = seq->name.l;
                        memcpy(s->name, seq->name.s, s->nlen);

                        this->readQueue_.push(s);
                        ksv = kseq_read(seq);
                      }

                    // destroy the parser and close the file
                    kseq_destroy(seq);
                    close(pfd.fd);
                }


                this->parsing_ = false;
            });
            return true;
        } else {
            return false;
        }

    }

bool StreamingReadParser::nextRead(ReadSeq*& seq) {
        while(parsing_ or !readQueue_.empty()) {
            if (readQueue_.try_pop(seq)) { return true; }
        }
        return false;
}

void StreamingReadParser::finishedWithRead(ReadSeq*& s) { seqContainerQueue_.push(s); }
