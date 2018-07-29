#include "squares.h"
#include "gtest/gtest.h"

using namespace squares;

// compare with Table 1 from paper
TEST(squares_test, paper)
{
    constexpr double alpha = 0.001;
    EXPECT_NEAR(pvalue(15.5, 5), alpha, 0.00002);
    EXPECT_NEAR(pvalue(23.8, 50), alpha, 0.00002);
    // from MC have larger error
    EXPECT_NEAR(pvalue(25.6, 100), alpha, 0.00008);
}

TEST(squares_test, mathematica)
{
    constexpr unsigned n = 20;
    auto T = {2., 5., 10., 20., 50.};
    auto P = {0.8936721808595665, 0.42457437866154357, 0.06934906413527009,
              0.0014159488909252227, 9.575188641974819e-9};
    auto eps = {1e-15, 1e-15, 1e-15, 1e-15, 1e-15};
    for (auto t = T.begin(), p = P.begin(), e=eps.begin(); t != T.end(); ++t, ++p, ++e)
        EXPECT_NEAR(pvalue(*t, n), *p, *e);
}

TEST(squares_test, critical)
{
    constexpr unsigned N = 50;

    EXPECT_NEAR(pvalue(19.645, 2*N), 0.01, 3e-5);
    EXPECT_NEAR(pvalue(15.34, 2*N), 0.05, 1e-4);
}
