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

#include "partitions.h"

#include <cassert>
#include <algorithm>
#include <iostream>

using namespace std;

Partition::Partition(UInt_t n) :
    c(n+1, 0),
    y(n+1, 0),
    n(n),
    h(1)
{
    assert(n > 0);

    y[0] = -1;
    c.at(0) = 1;
    y.at(1) = n;
    c[1] = 1;
}

Partition::Partition(UInt_t n, UInt_t k) :
    // have at most k distinct parts + one buffer value at index 0
    c(k+1, 0),
    y(k+1, 0),
    // also the index of the largest part
    n(n),
    h(0)
{
    assert(n > 0);
    assert(k > 0);
    assert(n >= k);

    y[0] = -1;
    c[0] = +1;

    // the largest part of any partition of n into k parts
    const auto maxPart = n-k+1;

    // define initial partition

    // If there is only one distinct part then maxPart divides n.  But
    // part(6,3) = 5+1 = 4+2 = 3+3 shows that it doesn't have to be
    // the first partition, rather it is always the last in our order
    // because it cannot be mutated further. So if first and last
    // agree, there is only one partition. From a quick glance at the
    // partition table or the recurrence relation defining it T(n, 1)
    // = T(n, n) = 1, T(n, k) = 0 (k>n), T(n, k) = T(n-1, k-1) +
    // T(n-k, k), this reduces the choices to k=1,n-1,n. But for
    // k=n-1, the partition is n = (n-2)*1 + 1*2 so has two distinct
    // parts if n > 2. For n=2, k=n-1=1 is no special case, and for
    // n=1 all options collapse.
    if (k == 1 || k == n) {
        y[1] = maxPart;
        c[1] = n / maxPart;
        h = 1;
    } else {
        // initialize to the partition where the largest part has the
        // maximum value and everything else is 1. There are two
        // distinct parts
        y[1] = 1;
        c[1] = k-1;
        y[2] = maxPart;
        c[2] = 1;
        h = 2;
    }
}

/*!
 * Code up the update step inside the while loop of algorithm Z from
 *
 * A. Zoghbi: Algorithms for generating integer partitions, Ottawa (1993)
 *
 * to generate the next integer partition of `n` in multiplicity
 * representation `\sum_i c_i * y_i` such that `c_i` is the
 * multiplicity of part `y_i`.
 *
 */
Partition& Partition::operator++()
{
    // running index
    Int_t i = h-1;
    // calculated part
    Int_t k = c[h];
    // calculated remainder
    Int_t r = c[h] * y[h];
    // calculate the remainder
    while (y[h] - y[i] < 2) {
        k += c[i];
        r += c[i] * y[i];
        --i;
    }

    // update current part when it equals 1
    if (c[i] == 1) {
        if (i != 0) {
            r += c[i] * y[i];
            ++y[i];
        } else {
            i = 1;
            y[i] = 1;
        }
    }
    // update current part != 1
    else {
        --c[i];
        r += y[i];
        ++i;
        y[i] = y[i-1] + 1;
    }

    // calculate next parts based on remainder left from previous update
    c[i] = k;
    r -= c[i] * y[i];
    h = i+1;

    // update last modified part if it'k the remainder
    if (r == y[i]) {
        ++c[i];
        h = i;
    }
    // add new part with multiplicity 1
    else {
        y[h] = r;
        c[h] = 1;
    }

    return *this;
}

bool Partition::operator==(const Partition& other) const
{
    if (n != other.n || h != other.h)
        return false;
    for (size_t i = 1; i <= h; ++i) {
        if (c[i] != other.c[i] || y[i] != other.y[i])
            return false;
    }
    return true;
}

AbstractPartitionGenerator::AbstractPartitionGenerator(const Partition& p) :
    done(false),
    p(p)
{

}

AbstractPartitionGenerator& AbstractPartitionGenerator::operator++()
{
    // could lead to segfault if we advance in the sequence of partitions
    assert(!done);

    // check if we reached the final partition
    if (final_partition())
        done = true;
    else
        ++p;

    return *this;
}

KPartitionGenerator::KPartitionGenerator(UInt_t n, UInt_t k) :
    AbstractPartitionGenerator(Partition(n, k))
{
}


KPartitionGenerator& KPartitionGenerator::operator++()
{
    AbstractPartitionGenerator::operator++();
    return *this;
}

bool KPartitionGenerator::final_partition() const
{
    return p.parts()[p.distinct_parts()] - p.parts()[1] <= 1;
}

PartitionGenerator::PartitionGenerator(UInt_t n) :
    AbstractPartitionGenerator(Partition(n))
{
}

PartitionGenerator& PartitionGenerator::operator++()
{
    AbstractPartitionGenerator::operator++();
    return *this;
}

bool PartitionGenerator::final_partition() const
{
    return UInt_t(p.mult()[1]) == p.number();
}

std::ostream& operator<<(std::ostream& out, const Partition& p)
{
    bool init = false;
    const auto h = p.distinct_parts();
    const auto& c = p.mult();
    const auto& y = p.parts();
    for(size_t i=1; i <= h; ++i) {
        out << (init? " + " : "") << c[i] << "*" << y[i];
        init = true;
    }
    return out;
}

// Local Variables:
// compile-command:"g++ -std=c++11 -g -O2 partitions.cxx -o partitions && ./partitions"
// End:
