// Copyright 2017 Frederik Beaujean <Frederik.Beaujean@lmu.de>

#include "pvalue.h"
#include "partitions.h"

#include <gsl/gsl_cdf.h>

#include <algorithm>
#include <cassert>
#include <cmath>

using namespace std;

// use ldouble for higher precision when adding lots of numbers
using ldouble = long double;

namespace
{
static std::vector<ldouble> log_factorial;

unsigned CacheFactorials(unsigned N)
{
    if (N < log_factorial.size())
        // log factorials have already been cached up to N
        return log_factorial.size();

    // reserve memory
    log_factorial.reserve(N);

    // add log(0!) if not there
    if (log_factorial.empty())
        log_factorial.push_back(0L);

    // calculate new log factorials
    for (unsigned i = log_factorial.size(); i <= N; ++i)
        log_factorial.push_back(log_factorial.back() + log(ldouble(i)));

    return log_factorial.size();
}

std::vector<ldouble> CacheChi2(double Tobs, unsigned N)
{
    assert(N>0);
    // to ease addressing, pad with zero element that, if used, should spoil any calculation
    std::vector<ldouble> res(N+1);
    res[0] = std::numeric_limits<ldouble>::quiet_NaN();

    // C style
#pragma omp parallel for shared(Tobs, N, res)
    for (size_t i = 1; i <= N; ++i)
        res[i] = log(ldouble(gsl_cdf_chisq_P(Tobs, i)));
    return res;
}

} // namespace

double runs_cumulative(const double Tobs, const unsigned N)
{
    CacheFactorials(N);

    // pretabulate chi2 cumulative: given N, we need P(Tobs|i) for i=1...N
    auto log_cumulative = CacheChi2(Tobs, N);

    // work on log scale to avoid overflows of Pochhammer symbol,
    // factorial, and the exponential. Use natural log
    ldouble poch = 0;
    // log(2^N-1): bit shift if N small enough, else neglect -1
    const ldouble logpow2N1 = (N<=63)? log((1ul << N) - 1) : N*log(2);

    // the p value
    ldouble p = 0;

// in tests for N=10 dynamic was better than schedule(static,2). Since
// part(r, M) is really different, each iteration can vary in time
// very much. Hyperthreading seemed to help a lot on my Intel K4770.
#pragma omp parallel for schedule(dynamic) private(poch) shared(log_cumulative, log_factorial) reduction(+:p)
    for (auto r=1ul; r <= N; ++r) {
        // const unsigned long Mmax = (r <= N-r+1)? r : N-r+1;
        const auto Mmax = min(r, N-r+1);
        poch = 0;
        for (auto M=1ul; M <= Mmax; ++M) {
          // compute Pochhammer iteratively
            // (...)_{M+1}/ (...)_M = N-r+1-M but to start the iteration we have to add 1
            poch += log(ldouble(N-r+2-M));

            // this factor is independent of the actual partition,
            // only depends on M,r,N
            const ldouble scale = poch - logpow2N1;

            // maintain sum over partitions
            ldouble ppi = 0;

            // visit all partitions, save ref to partition
            KPartitionGenerator g(r, M);
            auto& n = g->mult();
            auto& y = g->parts();

            for (; g; ++g) {
                const auto h = g->distinct_parts();

                // perform sum on log scale within a partition
                ldouble ppartition = 0;

                for (size_t l = 1; l <= h; ++l) {
                    ppartition += n[l] * log_cumulative[y[l]] - ::log_factorial[n[l]];
                }
                ppi += exp(ppartition);
            }

            // have to stay on linear scale
            p += exp(scale + log(ppi));
        }
    }
    assert(p < 1);

    return p;
}

double runs_pvalue(const double Tobs, const unsigned N)
{
    return 1 - runs_cumulative(Tobs, N);
}
