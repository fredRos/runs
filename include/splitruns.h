// Copyright (c) 2017 Frederik Beaujean (Frederik.Beaujean@lmu.de)

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

constexpr double splitruns_epsrel = 1e-10;
constexpr double splitruns_epsabs = 1e-15;

/**
 * F(Tobs | n*N). Evaluates F(Tobs | N) exactly, then approximates.
 */
double runs_split_cumulative(const double Tobs, const unsigned N, const double n, double epsrel=splitruns_epsrel, double epsabs=splitruns_epsabs);
double runs_split_pvalue(const double Tobs, const unsigned N, const double n, double epsrel=splitruns_epsrel, double epsabs=splitruns_epsabs);

/**
 * Compute \Delta correction term.
 */
double Delta(const double Tobs, const unsigned Nl, const unsigned Nr, double epsrel=splitruns_epsrel, double epsabs=splitruns_epsabs);

double h(const double chisq, const unsigned N);
double H(const double a, const double b, const unsigned N);
