#ifndef FRAGMENT_LIST_HPP
#define FRAGMENT_LIST_HPP

extern "C" {
#include "container.h"
}

class FragmentList {
    public:
        FragmentList();

         // Free all of the structures allocated for this fragment list
        ~FragmentList();

        void freeFragList(Container* frags);

        void computeBestChain_(Container* frags, double& maxScore);

        void computeBestChain();

       void addFragMatch(uint32_t refStart, uint32_t queryStart, uint32_t len);

        void addFragMatch(uint32_t refStart, uint32_t refEnd,
                         uint32_t queryStart, uint32_t queryEnd);

        void addFragMatchRC(uint32_t refStart, uint32_t queryStart, uint32_t len);

        void addFragMatchRC(uint32_t refStart, uint32_t refEnd,
                         uint32_t queryStart, uint32_t queryEnd);

        Container* fragments;
        Container* fragmentsRC;
        size_t numFrag;
        double bestHitScore;
};

#endif // FRAGMENT_LIST_HPP
