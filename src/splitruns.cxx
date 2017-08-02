#include "splitruns.h"
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
    double res = 0;
    double weight = 0.5;
    for (auto i = 1u; i <= N; ++i) {
        if (i < N)
            weight *= 0.5;
        res += weight * gsl_ran_chisq_pdf(chisq, i);
    }

    return res;
}

double H(const double a, const double b, const unsigned N)
{
    double res = 0;
    double weight = 0.5;
    for (auto i = 1u; i <= N; ++i) {
        if (i < N)
            weight *= 0.5;
        res += weight * (gsl_cdf_chisq_P(b, i) - gsl_cdf_chisq_P(a, i));
    }

    return res;
}

double integrand(double x, void* params)
{
    const IntegrandData & d = *static_cast<IntegrandData*>(params);

    return h(x, d.N) * H(d.Tobs - x, d.Tobs, d.N);
}

double Delta(const double Tobs, const unsigned N, double epsrel, double epsabs)
{
    // gsl numerical integration
    constexpr size_t limit = 1000;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc(limit);

    double result, error;

    gsl_function F;
    F.function = &integrand;
    ::IntegrandData data{Tobs, N};
    F.params = &data;

    gsl_integration_qag(&F, 0, Tobs, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, &result, &error);

    // printf ("result          = % .17e\n", result);
    // printf ("estimated error = % .6e\n", error);
    // printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    return result;
}

double runs_split_cumulative(const double Tobs, const unsigned N, const double n, double epsrel, double epsabs)
{
    const auto F = runs_cumulative(Tobs, N);
    const auto Fn1 = pow(F, n-1);
    return F * Fn1 / (1 + Delta(Tobs, N, epsrel, epsabs) * (1-Fn1) / (1-F));
}

double runs_split_pvalue(const double Tobs, const unsigned N, const double n, double epsrel, double epsabs)
{
    return 1-runs_split_cumulative(Tobs, N, n, epsrel, epsabs);
}
