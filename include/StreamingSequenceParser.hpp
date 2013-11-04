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
