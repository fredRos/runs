#include "splitruns.h"
#include "pvalue.h"
#include "gtest/gtest.h"
#include <gsl/gsl_cdf.h>
#include <cmath>

TEST(splitruns, split)
{
    constexpr double Tobs = 15.5;
    constexpr unsigned N = 12;
    constexpr unsigned n = 2;

    const double F = runs_cumulative(Tobs, N);
    const double approx = F*F / (1.0 + Delta(Tobs, N, N));

    EXPECT_NEAR(approx, runs_split_cumulative(Tobs, N, n), 1e-15);
}

void test_on_grid(const unsigned K, const unsigned N, const unsigned n, const double eps=2e-7)
{
    // check for absolute difference, same as relative error when all values near one.
    // agreement better when cumulative closer to one, or Tobs larger
    for (auto i = 1u; i < K; ++i) {
        const double Tobs = 20 + 2*i;
        EXPECT_NEAR(runs_split_cumulative(Tobs, N, n), runs_cumulative(Tobs, n*N), eps)
            << " at Tobs = " << Tobs << ", N = " << N << ", and n = " << n;
    }
}

TEST(splitruns, approx)
{
    // compare with full results
    constexpr unsigned K = 15;
    constexpr unsigned N = 40;
    test_on_grid(K, N, 2, 2e-7);
    // test_on_grid(K, N, 3, 4e-7);
}

TEST(splitruns, hH)
{
    constexpr double Tobs = 15.5;
    constexpr double x = 3.3;
    constexpr unsigned N = 12;

    // compare to mathematica
    EXPECT_NEAR(h(Tobs, N), 0.000373964, 1e-8);
    EXPECT_NEAR(H(Tobs-x, Tobs, N), 0.00245352, 1e-8);
    EXPECT_NEAR(Delta(Tobs, N, N), 0.00175994, 1e-8);
}

TEST(splitruns, cdf)
{
    EXPECT_NEAR(gsl_cdf_chisq_P(15.5, 12), 0.784775, 1e-6);
}

TEST(splitruns, long)
{
    // just use gtest to measure the time
    constexpr double Tobs = 32;
    constexpr unsigned N = 60;
    // constexpr unsigned n = 5;
    // runs_split_cumulative(Tobs, N, n, 1e-13, 1e-20);
    Delta(Tobs, N, N, 1e-4, 1e-20);
}

TEST(splitruns, interpolate)
{
    constexpr double Tobs = 32;

    std::cout << std::setprecision(15);
    double F71 = runs_split_cumulative(Tobs, 71, 5);
    double F88 = runs_split_cumulative(Tobs, 88, 4);
    double F89 = runs_split_cumulative(Tobs, 89, 4);
    double F100 = runs_split_cumulative(Tobs, 100, 3.55);
    std::cout << "F(32|5*71) = " << F71 << std::endl;
    std::cout << "1/3 F(32|4*88) + 3/4 F(32|4*89) = " << (F88 + 3*F89)/4 << std::endl;
    std::cout << "F(32|3.55*100) = " << F100 << std::endl;
    EXPECT_NEAR(F71, F100, 3e-9);
    EXPECT_NEAR(F71,(F88 + 3*F89)/4, 2e-9);
}

TEST(splitruns, bound_error)
{
    // the true correction is bounded by scaling by F(T) and F(2T)
    constexpr double Tobs = 8;
    EXPECT_NEAR(runs_cumulative(Tobs, 100), std::pow(runs_cumulative(Tobs, 50), 2) - runs_cumulative(1*Tobs, 100) * Delta(Tobs, 100, 100), 1e-3);
}

TEST(splitruns, 2dcorrection)
{
    constexpr double Tobs = 15.5;
    constexpr unsigned N = 20;
    const auto corr = full_correction(Tobs, N, N, 1e-7, 0.0, 20);

    // compare to mathematica when cumulative is ignored
    // EXPECT_NEAR(corr, 0.00175994, 1e-6);

    // the true result should lie within the bounds
    const auto delta = Delta(Tobs, N, N);
    const auto lo = runs_cumulative(Tobs, 2*N) * delta;
    const auto hi = runs_cumulative(2*Tobs, 2*N) * delta;

    std::cout << "lo = " << lo << ", full corr  = " << corr << ", hi = " << hi  << std::endl;

    EXPECT_LE(lo, corr);
    EXPECT_LE(corr, hi);
}

TEST(splitruns, paper_timing)
{
    constexpr double Tobs = 15.8;
    constexpr unsigned N = 96;
    // don't care about result, just want to measure the time
    runs_pvalue(Tobs, N);
}
