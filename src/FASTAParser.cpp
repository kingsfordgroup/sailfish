#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include "FASTAParser.hpp"
#include "Transcript.hpp"
#include "SailfishStringUtils.hpp"

extern "C" {
#include "kseq.h"
}

KSEQ_INIT(int, read)

FASTAParser::FASTAParser(const std::string& fname): fname_(fname) {}

void FASTAParser::populateTargets(std::vector<Transcript>& refs) {
    using std::string;
    using std::unordered_map;

    unordered_map<string, size_t> nameToID;
    for (size_t idx = 0; idx < refs.size(); ++idx) { nameToID[refs[idx].RefName] = idx; }
    
    FILE* fp = fopen(fname_.c_str(), "r");
    kseq_t* seq = kseq_init(fileno(fp));

    while (kseq_read(seq) >= 0) {
        std::string name(seq->name.s);
        auto it = nameToID.find(name);
        if (it == nameToID.end()) {
            std::cerr << "WARNING: Transcript " << name << " appears in the reference but did not appear in the BAM\n";
        } else {
            refs[it->second].Sequence = sailfish::stringtools::encodeSequenceInSAM(seq->seq.s, seq->seq.l);
        }
    }

    fclose(fp);
}

