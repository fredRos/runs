#include "squares_approx.h"
#include "squares.h"

#include "cubature.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_spline.h>

#include <cassert>
#include <cmath>
#include <stdexcept>

namespace {
    struct IntegrandData
    {
        double Tobs;
        unsigned Nl, Nr;
    };

    struct CubaIntegrandData
    {
        double Tobs;
        unsigned Nl, Nr;
        unsigned counter;
        gsl_spline* spline;
        gsl_interp_accel* acc;
    };
}

namespace squares
{

double h(const double chisq, const unsigned N)
{
    double res = 0;
    double weight = 0.5;
    for (auto i = 1u; i <= N; ++i)
    {
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
    for (auto i = 1u; i <= N; ++i)
    {
        if (i < N)
            weight *= 0.5;
        res += weight * (gsl_cdf_chisq_P(b, i) - gsl_cdf_chisq_P(a, i));
    }

    return res;
}

double gsl_integrand(double x, void *params)
{
    const IntegrandData &d = *static_cast<IntegrandData *>(params);

    return h(x, d.Nl) * H(d.Tobs - x, d.Tobs, d.Nr);
}

double Delta(const double Tobs, const unsigned Nl, const unsigned Nr, double epsrel, double epsabs)
{
    // gsl numerical integration
    constexpr size_t limit = 1000;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);

    double result, error;

    gsl_function F;
    F.function = &gsl_integrand;
    ::IntegrandData data{Tobs, Nl, Nr};
    F.params = &data;

    gsl_integration_qag(&F, 0, Tobs, epsabs, epsrel, limit, GSL_INTEG_GAUSS21, w, &result, &error);

    // printf ("result          = % .17e\n", result);
    // printf ("estimated error = % .6e\n", error);
    // printf ("intervals       = %zu\n", w->size);

    gsl_integration_workspace_free(w);

    return result;
}

int cubature_integrand(unsigned /*ndim*/, const double *uv, void *data, unsigned /*fdim*/, double *fval)
{
    CubaIntegrandData &d = *static_cast<CubaIntegrandData *>(data);
    ++d.counter;

    /*  transform from unit square to triangle with vertices (Tobs, Tobs), (Tobs, 0), (0, Tobs)  */
    const auto u = uv[0];
    const auto v = uv[1];
    const auto x = d.Tobs * u * (1 - v) + d.Tobs * (1 - u);
    const auto y = d.Tobs * u * v + d.Tobs * (1 - u);
    const auto jac = d.Tobs * d.Tobs * u;

    *fval = jac * h(x, d.Nl) * h(y, d.Nr);
    *fval *= (d.spline && d.acc) ? gsl_spline_eval(d.spline, x + y, d.acc) : cumulative(d.Tobs, d.Nl + d.Nr);

    return 0;
}

double full_correction(const double Tobs,
                       const unsigned Nl,
                       const unsigned Nr,
                       double epsrel,
                       double epsabs,
                       unsigned ninterp)
{
    ::CubaIntegrandData data{Tobs, Nl, Nr, 0, nullptr, nullptr};

    /* evaluate F on a grid (Tobs, ..., 2*Tobs) for interpolation */
    gsl_interp_accel *acc = nullptr;
    gsl_spline *spline = nullptr;
    double *x = nullptr;
    double *y = nullptr;

    if (ninterp >= 2)
    {
        acc = gsl_interp_accel_alloc();
        // gsl_interp_steffen would be differentiable at grid points but is only available in newer versions of the GSL. We want monotonicity because if F is, too
        spline = gsl_spline_alloc(gsl_interp_linear, ninterp);

        x = new double[ninterp];
        y = new double[ninterp];

        for (auto i = 0u; i < ninterp; ++i)
        {
            x[i] = Tobs + i * Tobs / (ninterp - 1);
            y[i] = cumulative(x[i], Nl + Nr);
        }

        gsl_spline_init(spline, x, y, ninterp);
    }
    data.acc = acc;
    data.spline = spline;

    constexpr unsigned fdim = 1;
    constexpr unsigned dim = 2;
    constexpr double uvmin[dim] = {0.0, 0.0};
    constexpr double uvmax[dim] = {1.0, 1.0};
    constexpr size_t maxEval = 10000;
    double res;
    double err;

    if (hcubature(fdim, cubature_integrand, &data, dim, uvmin, uvmax,
                  maxEval, epsabs, epsrel, ERROR_L2, &res, &err))
    {
        throw std::runtime_error("hcubature failed");
    }
    // printf("Computed integral = %0.10g +/- %g with %d calls\n", res, err, data.counter);

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    delete x;
    delete y;

    return res;
}

double approx_cumulative(const double Tobs, const unsigned N, const double n, double epsrel, double epsabs)
{
    const auto F = cumulative(Tobs, N);
    const auto Fn1 = pow(F / (1 + Delta(Tobs, N, N, epsrel, epsabs)), n - 1);
    return F * Fn1;
}

double approx_pvalue(const double Tobs, const unsigned N, const double n, double epsrel, double epsabs)
{
    return 1 - approx_cumulative(Tobs, N, n, epsrel, epsabs);
}

}
