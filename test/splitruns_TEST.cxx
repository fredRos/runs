#include "splitruns.h"
#include "pvalue.h"
#include "gtest/gtest.h"
#include <gsl/gsl_cdf.h>

TEST(splitruns, split)
{
    constexpr double Tobs = 15.5;
    constexpr unsigned N = 12;
    constexpr unsigned n = 2;

    const double F = runs_cumulative(Tobs, N);
    const double approx = F*F / (1.0 + Delta(Tobs, N));

    EXPECT_NEAR(approx, runs_split_cumulative(Tobs, N, n), 1e-16);
}

TEST(splitruns, approx)
{
    // compare with full results
    constexpr unsigned N = 25;
    constexpr unsigned n = 2;

    // check for absolute difference here but relative error would be
    // more interesting. Not implemented in gtest
    for (auto i = 1; i < 20; ++i) {
        const double Tobs = 20 + 2*i;
        EXPECT_NEAR(runs_split_cumulative(Tobs, N, n), runs_cumulative(Tobs, n*N), 2e-7)
            << " at Tobs = " << Tobs;
    }
}

TEST(splitruns, hH)
{
    constexpr double Tobs = 15.5;
    constexpr double x = 3.3;
    constexpr unsigned N = 12;

    // compare to mathematica
    EXPECT_NEAR(h(Tobs, N), 0.000373964, 1e-8);
    EXPECT_NEAR(H(Tobs-x, Tobs, N), 0.00245352, 1e-8);
    EXPECT_NEAR(Delta(Tobs, N), 0.00175994, 1e-8);
}

TEST(splitruns, cdf)
{
    EXPECT_NEAR(gsl_cdf_chisq_P(15.5, 12), 0.784775, 1e-6);
}
