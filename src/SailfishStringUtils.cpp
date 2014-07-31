#include "SailfishStringUtils.hpp"
#include <cstdint>

uint8_t* sailfish::stringtools::encodeSequenceInSAM(const char* src, size_t len) {
    uint8_t* target = new uint8_t[len];
    for(size_t i = 0; i < len; ++i) {
        size_t byte = i >> 1;
        size_t nibble = i & 0x1;
        if (nibble) {
            target[byte] |= charToSamEncode[src[i]];
        } else {
            target[byte] |= (charToSamEncode[src[i]] << 4);
        }
    }
    return target;
}

