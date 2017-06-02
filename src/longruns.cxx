#include "longruns.h"
#include "pvalue.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>

#include <cassert>
#include <cmath>

namespace {
    struct IntegrandData
    {
        double Tobs;
        unsigned N;
    };
}

double h(const double chisq, const unsigned N)
{
    // to avoid overflow
    assert(N < 64);

    double res = 0;
    for (auto i = 1u; i <= N; ++i) {
        unsigned K = (i == N) ? N : i+1;
        res += 1.0 / (1 << K) * gsl_ran_chisq_pdf(chisq, i);
    }

    return res;
}

double H(const double a, const double b, const unsigned N)
{
    // to avoid overflow
    assert(N < 64);

    double res = 0;
    for (auto i = 1u; i <= N; ++i) {
        // power of 2
        unsigned K = (i == N) ? N : i+1;
        double weight = 1.0 / (1ul << K);
        res += weight * (gsl_cdf_chisq_P(b, i) - gsl_cdf_chisq_P(a, i));
    }

    return res;
}

double integrand(double x, void* params)
{
    const IntegrandData & d = *static_cast<IntegrandData*>(params);

    return h(x, d.N) * H(d.Tobs - x, d.Tobs, d.N);
}

double Delta(const double Tobs, const unsigned N)
{
    // gsl numerical integration
    constexpr size_t limit = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

    double result, error;

    gsl_function F;
    F.function = &integrand;
    ::IntegrandData data{Tobs, N};
    F.params = &data;

    gsl_integration_qags (&F, 0, Tobs, 0, 1e-8, limit,
                          w, &result, &error);

    // printf ("result          = % .6e\n", result);
    // printf ("estimated error = % .6e\n", error);
    // printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    return result;
}

double runs_split_cumulative(const double Tobs, const unsigned N, const unsigned n)
{
    const double F = runs_cumulative(Tobs, N);
    return pow(F, n) / (1 + (n-1) * Delta(Tobs, N));
}

double runs_split_pvalue(const double Tobs, const unsigned N, const unsigned n)
{
    return 1-runs_split_cumulative(Tobs, N, n);
}
