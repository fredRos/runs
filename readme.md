Runs
====

This project contains a mathematica and a C++ implementation of the
distribution of the runs statistic defined in

Frederik Beaujean and Allen Caldwell. “A Test Statistic for Weighted
Runs.” Journal of Statistical Planning and Inference 141, no. 11
(November 2011):
3437–46. [doi:10.1016/j.jspi.2011.04.022](http://dx.doi.org/10.1016/j.jspi.2011.04.022) [arXiv:1005.3233](http://arxiv.org/abs/1005.3233).

mathematica
------------

The reference implementation is in the `mathematica` package
`RunsWeightedPackage`. To use it from a `mathematica` notebook,
download the package to a directory `/package/dir` on your computer and do

    SetDirectory["/package/dir"]
    Needs["RunsWeightedPackage`"]

All exported commands in the package are accessible via

    ?"RunsWeightedPackage`*"

The most important commands to compute the p value for the runs
statistic are

    # compute the p value using the exact expression that involves
    # summing over integer partitions
    n = 10
    pValueRuns[3.3, n]

    # use Monte Carlo experiments to approximate the p value. Much
    # faster for n > 80 than the exact expression
    n = 500
    pValueRunsMC[35.3, n]

    # calculate the runs statistic from a list of values interpreted
    # as independent samples from a standard normal distribution
    runsSuccess
    runsSuccess[{-1,1,3,-2}]==10

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

    #include <partition.h>
    #include <iostream>

    # print all partitions of 6 into any number of parts
    for (PartitionGenerator gen(6); gen; ++gen)
        std::cout << *gen << std::endl;

    # print all partitions of 6 into 3 parts
    for (KPartitionGenerator gen(6, 3); gen; ++gen)
        std::cout << *gen << std::endl;

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

### runs statistics

`Tobs` denotes the value of the runs test statistic; i.e., the largest
`\chi^2` of any run of consecutive successes (above expectation) in a
sequence of `N` independent trials with Gaussian uncertainty. Then the
cumulative distribution `P(T < Tobs | N)` and the p value `P(T >= Tobs
| N)` are available as

    include "pvalue.h"

    runs_cumulative(Tobs, N);
    runs_pvalue(Tobs, N);

`openMP` helps as the speed-up of evaluating `runs_cumulative` for
large `N>50` scales linearly with the number of physical cores and
even benefits from hyperthreading. We tested on Intel Core i7-4600U
with gcc4.9 and Intel Core i7-4770 with gcc 5.4.

build instructions
------------------

`cmake`, `c++11`, and `GSL` are required. `OpenMP` is optional but recommended.

`google test` is needed for the tests and downloaded automatically when
building the first time.

    git clone https://github.com/fredRos/runs.git
    cd runs
    mkdir build
    cd build
    cmake ..
    make
    OMP_NUM_THREADS=4 ./runs_test

citing
------

If you use this code in an academic setting, please cite this
reference.

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