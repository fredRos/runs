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

#include <ostream>
#include <iterator>
#include <vector>

using UInt_t = unsigned;
using Int_t = int;
using vec = std::vector<Int_t>;

/*!
 * Represent a partition of a positive number in multiplicity representation.
*/
class Partition
{
public:
    Partition(UInt_t n);
    Partition(UInt_t n, UInt_t k);
    Partition& operator++();
    bool operator==(const Partition&) const;
    const vec& mult() const { return c; }
    const vec& parts() const { return y; }
    UInt_t distinct_parts() const { return h; }
    UInt_t number() const { return n; }

private:
    vec c;     /// multiplicity
    vec y;     /// part
    const UInt_t n; /// Partition of n
    UInt_t h;  /// number of distinct parts
};

/*!
 * Generate all integer partitions in multiplicity representation.
 *
 * Based on the pseudocode of Algorithm Z in A. Zoghbi: Algorithms for
 * generating integer partitions, Ottawa (1993), which generates all
 * partitions of `n`. I made the necessary modifications to partition
 * into `k` parts.
 *
 * Example:
 * for (PartitionGenerator gen(6); gen; ++gen)
 *       std::cout << *gen << std::endl;
 *
 * Output:
 * 1*6
 * 1*1 + 1*5
 * 1*2 + 1*4
 * 2*3
 * 2*1 + 1*4
 * 1*1 + 1*2 + 1*3
 * 3*2
 * 3*1 + 1*3
 * 2*1 + 2*2
 * 4*1 + 1*2
 * 6*1
 */
class AbstractPartitionGenerator : public std::iterator<std::input_iterator_tag, Partition>
{
public:
    explicit operator bool() const { return !done; }

    // We are happy with most default types for std::iterator but we
    // don't want the caller to modify a partition
    using reference = value_type const&;
    using pointer = value_type const*;

    reference operator*() const { return p; }
    pointer operator->() const { return &p; }

    /// Increment to next partition. Fails if already done; i.e. bool(this) == false
    AbstractPartitionGenerator& operator++();

protected:
    AbstractPartitionGenerator(const Partition&);
    virtual bool final_partition() const = 0;
    bool done;
    value_type p;
};

/// Generate all partitions of `n` into `k` parts.
class KPartitionGenerator : public virtual AbstractPartitionGenerator
{
public:
    KPartitionGenerator(UInt_t n, UInt_t k);

    KPartitionGenerator& operator++();
protected:
    virtual bool final_partition() const override;
};

/// Generate all partitions of `n`.
class PartitionGenerator : public virtual AbstractPartitionGenerator
{
public:
    PartitionGenerator(UInt_t n);

    PartitionGenerator& operator++();
protected:
    virtual bool final_partition() const override;
};


std::ostream& operator<<(std::ostream&, const Partition&);
