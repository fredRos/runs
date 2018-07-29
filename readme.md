[![Build Status](https://travis-ci.org/fredRos/runs.svg?branch=master)](https://travis-ci.org/fredRos/runs) [![DOI](https://zenodo.org/badge/87222105.svg)](https://zenodo.org/badge/latestdoi/87222105)

Weighted-runs (SQUARES) statistic
===========================================

This project contains a mathematica and a C++ implementation of the cumulative
distribution function of the weighted-runs statistic originally
defined in

Frederik Beaujean and Allen Caldwell. _A Test Statistic for Weighted
Runs._ Journal of Statistical Planning and Inference 141, no. 11
(November 2011):
3437–46. [doi:10.1016/j.jspi.2011.04.022](http://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](http://arxiv.org/abs/1005.3233)

We derived an approximation to be able to compute the cumulative also for large number of observations in

Frederik Beaujean and Allen Caldwell. _Is the bump significant? An axion-search example_ [arXiv:1710.06642](http://arxiv.org/abs/1710.06642)

where we renamed the weighted-runs statistic to the SQUARES statistic.

mathematica
------------

The reference implementation is in the `mathematica` package
`RunsWeightedPackage`. To use it from a `mathematica` notebook,
download the package to a directory `/package/dir` on your computer and do

``` mathematica
SetDirectory["/package/dir"]
Needs["RunsWeightedPackage`"]
```

All exported commands in the package are accessible via

``` mathematica
?"RunsWeightedPackage`*"
```

The most important commands to compute the p value for the runs
statistic are

``` mathematica
(* compute the p value using the exact expression that involves summing over integer partitions *)
n = 10
pValueRuns[3.3, n]

(* use Monte Carlo experiments to approximate the p value. Much faster for n > 80 than the exact expression *)
n = 500
pValueRunsMC[35.3, n]

(* calculate the runs statistic from a list of values interpreted as independent samples from a standard normal distribution *)
runsSuccess
runsSuccess[{-1,1,3,-2}]==10
```

c++
---

### partitions

The code features a standalone implementation of generating all
integer partitions of `n` in the multiplicity representation based on
the pseudocode of Algorithm Z in

A. Zoghbi: Algorithms for generating integer partitions, Ottawa
(1993), http://www.ruor.uottawa.ca/handle/10393/6506.

I made the necessary modifications to partition `n` into exactly `k`
parts as well. Some simple examples

``` c++
#include <partitions.h>
#include <iostream>

using namespace partitions;

// print all partitions of 6 into any number of parts
for (PartitionGenerator gen(6); gen; ++gen)
    std::cout << *gen << std::endl;

// print all partitions of 6 into 3 parts
for (KPartitionGenerator gen(6, 3); gen; ++gen)
    std::cout << *gen << std::endl;
```

`gen` is an iterator to a `Partition` and permits visiting all
partitions while keeping memory allocations to a minimum. A
`Partition` is essentially a triple `c, y, h` (in Zoghbi's notation) where

    n = \sum_{i=1}^h c_i * y_i

* `c = Partition::mult()` refers to a vector with the multiplicities,
* `y = Partition::parts()` refers to a vector with distinct parts,
* `h = Partition::distinct_parts()` refers to the number of distinct parts.

To reduce the memory allocations, the generator updates the partition
in place and the vectors `c,y` are long enough to hold the partition
with the maximum number of parts and a buffer element. For example,
the first partition of 6 into 3 parts is

    6 = 4+1+1
      = 2*1 + 1*4

In the code, this becomes

    h == 2

    c[1] == 2
    y[1] == 1

    c[2] == 1
    y[2] == 4

Ignore the first element `c[0], y[0]` and do not read beyond
`c[h],y[h]`!

### SQUARES statistic

`Tobs` denotes the value of the SQUARES test statistic; i.e., the largest
`\chi^2` of any run of consecutive successes (above expectation) in a
sequence of `N` independent trials with Gaussian uncertainty. Then the
cumulative distribution `P(T < Tobs | N)` and the p value `P(T >= Tobs| N)` are available as

``` c++
include "squares.h"

squares::cumulative(Tobs, N);
squares::pvalue(Tobs, N);
```

`openMP` helps as the speed-up of evaluating `squares::cumulative` for
large `N>50` scales linearly with the number of physical cores and
even benefits from hyperthreading. 

### split runs

For large `N`, the number of terms in the exact expressions scales like
`exp(N^1/2)/N` and quickly grows too large. We implement an approximate formula for
`n*N`, where for example `N = 100` is computed exactly and `n` may be 1 or >>
`N` and need not even be an integer.

```c++
#include "squares_approx.h"

squares::approx_cumulative(Tobs, N, n);
squares::approx_pvalue(Tobs, N, n);
```

The approximation involves a 1D numerical integration whose relative and absolute
target precision can be set as additional arguments, for example

```c++
approx_cumulative(Tobs, N, n, epsrel, epsabs)
```

For the exact meaning of these parameters consult
the
[GSL manual](https://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration-Introduction.html).

build instructions
------------------

`cmake`, `c++11`, and `GSL` are required. `OpenMP` is optional but recommended.

`google test` is needed for the tests and downloaded automatically when
building the first time.

    git clone https://github.com/fredRos/runs.git
    cd runs
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make
    OMP_NUM_THREADS=4 ./runs_test

To specify where the library and headers should be installed, use

    cmake -DCMAKE_INSTALL_PREFIX=/tmp/runs ..

To install

    make install

To link to this in your own code, just link to the library and include
the headers from whereever you chose to install them. For example if
you installed to `/tmp/runs/`, compile, link, and run your code in `runstest.cxx`

``` c++
#include "squares.h"
#include <iostream>

int main()
{
    double Tobs = 3.3;
    unsigned N = 10;
    std::cout << "F(Tobs = " << Tobs << " | N = " << N << ") = "
              << squares::cumulative(Tobs, N) << std::endl;

    return 0;
}
```

like this

    export LD_LIBRARY_PATH=/tmp/runs/lib/

    # with openmp: gcc >= 4.2, clang >= 3.8
    g++ -fopenmp runstest.cxx -I/tmp/runs/include -L/tmp/runs/lib -lruns -lgsl -lblas && ./a.out

    # without openmp: usually on a mac
    g++ runstest.cxx -I/tmp/runs/include -L/tmp/runs/lib -lruns -lgsl -lblas && ./a.out

testing
-------

In the `build` directory, call

    ./runs_test

To select individual tests

    ./runs_test --gtest_filter='splitruns.*'

Help on available options

    build/runs_test -h

strong scaling
-------

The C++ implementation of `squares::cumulative` and thus `squares::pvalue` benefits from
hyperthreading. On an Intel Core i7-4770 with four cores and a maximum frequency
of 3.4 GHz with gcc 5.4 and release mode, we observed the following run
times for the unit test `OMP_NUM_THREADS=n ./runs_test
--gtest_filter=splitruns.paper_timing` that computes `squares::pvalue(T, N=96)`

|     n |     time / ms |     speed up |
| :---: | :-----------: | :----------: |
|     1 |          9400 |            1 |
|     2 |          4750 |            2 |
|     4 |          2450 |          3.8 |
|     8 |          1750 |          5.4 |

This is a nice example where hyperthreading brings a noticeable improvement
beyond the number of physical cores.

citing
------

If you use this code in an academic setting, please cite these 
references.

    @article{beaujean2011test,
    title={A test statistic for weighted runs},
    author={Beaujean, Frederik and Caldwell, Allen},
    journal={Journal of Statistical Planning and Inference},
    volume={141},
    number={11},
    pages={3437--3446},
    year={2011},
    publisher={Elsevier}
    }
    @article{Beaujean:2017eyq,
      author         = "Beaujean, Frederik and Caldwell, Allen and Reimann, Olaf",
      title          = "{Is the bump significant? An axion-search example}",
      year           = "2017",
      eprint         = "1710.06642",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ex",
      SLACcitation   = "%%CITATION = ARXIV:1710.06642;%%"
}

license
-------

The code is released under the MIT license, see the `LICENSE` file. It comes
bundled with parts of the [cubature](https://github.com/stevengj/cubature)
package that is under the GPLv3.
