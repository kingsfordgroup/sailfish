/**
 *  This implementation for keeping an empirical distribution and
 *  querying the PDF and CDF is taken (and ever-so-slightly-modified) from
 *  https://raw.githubusercontent.com/dcjones/isolator/master/src/emp_dist.hpp.
 *  The original author is Daniel C. Jones; NOT me (Rob Patro).
 **/
#ifndef EMPIRICAL_DISTRIBUTION_HPP
#define EMPIRICAL_DISTRIBUTION_HPP

#include <climits>
#include <vector>


class EmpiricalDistribution
{
    public:
        EmpiricalDistribution(EmpiricalDistribution&);

        /* Construct a emperical distribution from n observations stored in xs.
         *
         * Input in run-length encoded samples, which must be in sorted order.
         *
         * Args:
         *   vals: An array of unique observations.
         *   lens: Nuber of occurances for each observation.
         *   n: Length of vals and lens.
         */
        EmpiricalDistribution(const std::vector<uint32_t>& vals,
                const std::vector<uint32_t>& count);

        /* Compute the median of the distribution. */
        float median() const;

        /* Probability density function. */
        float pdf(unsigned int x) const;

        /* Comulative p robabability. */
        float cdf(unsigned int x) const;

        /* The minimum observed value. */
        uint32_t minValue() const;

        /* The maximum observed value. */
        uint32_t maxValue() const;

        /* True if there are observations and false otherwise */
        bool valid() const;

    private:
        std::vector<float> pdfvals;
        std::vector<float> cdfvals;

        /* Precomputed median */
        float med;

        /* Min and Max Values */
        uint32_t minVal;
        uint32_t maxVal;
};



#endif
