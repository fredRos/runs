#include "longruns.h"
#include "pvalue.h"
#include "gtest/gtest.h"

TEST(longruns, paper)
{
    constexpr double Tobs = 15;
    constexpr unsigned N = 20;
    const double Delta = delta(Tobs, N);

    // compare with full results
    constexpr unsigned n = 2;
    const double target = runs_cumulative(Tobs, n*N);
    const double F = runs_cumulative(Tobs, N);
    const double approx = F*F / (1.0 + Delta);
    EXPECT_NEAR(approx, target, 1e-10);
}
