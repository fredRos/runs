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

    EXPECT_NEAR(approx, runs_split_cumulative(Tobs, N, n), 1e-15);
}

void test_on_grid(const unsigned K, const unsigned N, const unsigned n, const double eps=2e-7)
{
    // check for absolute difference here but relative error would be
    // more interesting. Not implemented in gtest
    for (auto i = 1u; i < K; ++i) {
        const double Tobs = 20 + 2*i;
        EXPECT_NEAR(runs_split_cumulative(Tobs, N, n), runs_cumulative(Tobs, n*N), eps)
            << " at Tobs = " << Tobs;
    }
}

TEST(splitruns, approx)
{
    // compare with full results
    constexpr unsigned K = 20;
    test_on_grid(K, 25, 2);
    test_on_grid(K, 25, 3, 4e-7);
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

TEST(splitruns, long)
{
    // just use gtest to measure the time
    constexpr double Tobs = 32;
    constexpr unsigned N = 60;
    // constexpr unsigned n = 5;
    // runs_split_cumulative(Tobs, N, n, 1e-13, 1e-20);
    Delta(Tobs, N, 1e-4, 1e-20);
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
    EXPECT_NEAR(F71, F100, 1e-9);
    EXPECT_NEAR(F71,(F88 + 3*F89)/4, 1e-9);
}
