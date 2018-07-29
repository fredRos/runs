// Copyright (c) 2018 Frederik Beaujean <beaujean@mpp.mpg.de>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

namespace squares
{

constexpr double EPSREL = 1e-10;
constexpr double EPSABS = 1e-15;

/**
 * F(Tobs | n*N). Evaluates F(Tobs | N) exactly, then approximates.
 */
double approx_cumulative(const double Tobs,
                         const unsigned N,
                         const double n,
                         double epsrel = EPSREL,
                         double epsabs = EPSABS);
double approx_pvalue(const double Tobs,
                     const unsigned N,
                     const double n,
                     double epsrel = EPSREL,
                     double epsabs = EPSABS);

/**
 * Compute \Delta correction term.
 */
double Delta(const double Tobs,
             const unsigned Nl,
             const unsigned Nr,
             double epsrel = EPSREL,
             double epsabs = EPSABS);

/**
 * Compute full correction w/o factoring out the cumulative.
 *
 * The 2D numerical integral is done by hcubature, which interprets
 * the relative and absolute precision `epsrel` and
 * `epsabs`. Expensive calls to `runs_cumulative` can be done once up
 * front and in the 2D integration, a linear interpolation is used in
 * place of `runs_cumulative` if the number of interpolation points
 * `ninterp >= 2`.
 */
double full_correction(const double Tobs,
                       const unsigned Nl,
                       const unsigned Nr,
                       double epsrel = EPSREL,
                       double epsabs = EPSABS,
                       unsigned ninterp = 0);

double h(const double chisq, const unsigned N);
double H(const double a, const double b, const unsigned N);

}
