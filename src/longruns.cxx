#include "longruns.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>

namespace {
    struct IntegrandData
    {
        double Tobs;
        unsigned N;
    };

    double h(const double chisq, const unsigned N)
    {
        double res = 0;
        // TODO could speed up a bit by multiplying previous term
        for (auto i = 1u; i < N; ++i) {
            res += 1.0 / (1 << (i+1)) * gsl_ran_chisq_pdf(chisq, i);
        }
        res += 1.0 / (1ul << N) * gsl_ran_chisq_pdf(chisq, N);

        return res;
    }

    double H(const double a, const double b, const unsigned N)
    {
        double res = 0;
        // TODO could speed up a bit by multiplying previous term
        for (auto i = 1u; i < N; ++i) {
            res += 1.0 / (1ul << (i+1)) * (gsl_cdf_chisq_P(b, i) - gsl_cdf_chisq_P(a, i));
        }
        res += 1.0 / (1 << N) * (gsl_cdf_chisq_P(b, N) - gsl_cdf_chisq_P(a, N));
        return res;
    }
}

double integrand(double x, void* params)
{
    const IntegrandData & d = *static_cast<IntegrandData*>(params);

    return h(x, d.N) * H(d.Tobs - x, d.Tobs, d.N);
}

double delta(const double Tobs, const unsigned N)
{
    // gsl numerical integration
    gsl_integration_workspace * w
        = gsl_integration_workspace_alloc (1000);

    double result, error;

    gsl_function F;
    F.function = &integrand;
    ::IntegrandData data{Tobs, N};
    F.params = &data;

    gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
                          w, &result, &error);

    printf ("result          = % .6e\n", result);
    printf ("estimated error = % .6e\n", error);
    printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free (w);

    return result;
}
