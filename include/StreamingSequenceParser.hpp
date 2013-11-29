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


#ifndef __STREAMING_READ_PARSER__
#define __STREAMING_READ_PARSER__

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread>
#include <atomic>
#include <iostream>
#include "unistd.h"
#include "fcntl.h"

extern "C" {
	#include "kseq.h"
}

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>
#include "tbb/concurrent_queue.h"

namespace bfs = boost::filesystem;

struct ReadSeq {
    char* seq = nullptr;
    size_t len = 0;
    char* name = nullptr;
    size_t nlen = 0;
};

class StreamingReadParser {
public:
    StreamingReadParser( std::vector<bfs::path>& files );
    ~StreamingReadParser();
    bool start();
    bool nextRead(ReadSeq*& seq);
    void finishedWithRead(ReadSeq*& s);

private:
    std::vector<bfs::path>& inputStreams_;
    bool parsing_;
    std::thread* parsingThread_;
    tbb::concurrent_bounded_queue<ReadSeq*> readQueue_, seqContainerQueue_;
    ReadSeq* readStructs_;
    const size_t queueCapacity_ = 2000000;
};

//#include "Parser.cpp"

#endif // __STREAMING_READ_PARSER__
