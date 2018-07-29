// Copyright 2018 Frederik Beaujean <beaujean@mpp.mpg.de>

#pragma once

namespace squares
{

/*!
 * Compute the cumulative distribution of the weighted-runs, or SQUARES, statistic `T`.
 *
 * The calculation implements Eqns. (16) and (17) from
 *
 * Frederik Beaujean and Allen Caldwell. “A Test Statistic for
 * Weighted Runs.” Journal of Statistical Planning and Inference 141,
 * no. 11 (November 2011): 3437–46. doi:10.1016/j.jspi.2011.04.022
 * http://arxiv.org/abs/1005.3233.
 *
 * @arg Tobs The value of the test statistic for the observed data set;
 * i.e., the largest \chi^2 of any run of consecutive observed values
 * above the expectation.
 * @arg N The total number of data points.
 */
double cumulative(const double Tobs, const unsigned N);
double pvalue(const double Tobs, const unsigned N);

}