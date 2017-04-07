#include "partitions.h"
#include "gtest/gtest.h"

#include <iostream>

template<class G>
void check(G& g, std::initializer_list<Int_t> c, std::initializer_list<Int_t> y)
{
    // check lengths of c and y match
    const auto h = c.size();
    ASSERT_EQ(h, y.size());

    // from 1 to h match values
    auto cit = c.begin();
    auto yit = y.begin();
    for (size_t i=1; i <= h; ++i, ++cit, ++yit) {
        EXPECT_EQ(g->mult().at(i), *cit) << " at index " << i;
        EXPECT_EQ(g->parts().at(i), *yit) << " at index " << i;
    }

    // now update partition to prepare next check
    ++g;
}

TEST(partitions_test, partition)
{
    // check partitions of n
    constexpr UInt_t n = 6;

    PartitionGenerator g(n);

    check(g, {1}, {6});
    check(g, {1,1}, {1,5});
    check(g, {1,1}, {2,4});
    check(g, {2}, {3});
    check(g, {2,1}, {1,4});
    check(g, {1,1,1}, {1,2,3});
    check(g, {3}, {2});
    check(g, {3,1}, {1,3});
    check(g, {2,2}, {1,2});
    check(g, {4,1}, {1,2});
    check(g, {6}, {1});
}

TEST(partitions_test, kpartition)
{
    KPartitionGenerator g(6, 3);

    check(g, {2,1}, {1,4});
    check(g, {1,1,1}, {1,2,3});
    check(g, {3}, {2});
}
