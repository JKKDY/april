#include <gtest/gtest.h>
#include <vector>
#include <ranges>
#include <algorithm>
#include <numeric>

#include "april/math/range.hpp"

using namespace april::math;

// ---------------------
// Constructors & Basics
// ---------------------

TEST(RangeTest, DefaultConstructor_IsEmpty) {
    constexpr Range r;
    static_assert(r.empty(), "Default constructed range must be empty in constexpr context");
    EXPECT_TRUE(r.empty());
    EXPECT_EQ(r.size(), 0u);
    EXPECT_EQ(r.start, 0u);
    EXPECT_EQ(r.stop, 0u);
}

TEST(RangeTest, StartStopConstructor_SetsBounds) {
    constexpr Range r(10, 20);
    EXPECT_EQ(r.start, 10u);
    EXPECT_EQ(r.stop, 20u);
    EXPECT_EQ(r.size(), 10u);
    EXPECT_FALSE(r.empty());
}

TEST(RangeTest, StartStopConstructor_ClampsStopToStart) {
    // If stop < start, size should be 0, start should be preserved, stop clamped to start
    constexpr Range r(10, 5);
    EXPECT_EQ(r.start, 10u);
    EXPECT_EQ(r.stop, 10u); // std::max(10, 5) -> 10
    EXPECT_TRUE(r.empty());
}

TEST(RangeTest, PairConstructor_WorksWithIntegers) {
    std::pair p{5, 15};
    Range r(p);
    EXPECT_EQ(r.start, 5u);
    EXPECT_EQ(r.stop, 15u);
    EXPECT_EQ(r.size(), 10u);
}



// -----------------
// Accessors & Logic
// -----------------

TEST(RangeTest, Contains_ReturnsTrueForValuesInRange) {
    Range r(10, 20);
    EXPECT_FALSE(r.contains(9));
    EXPECT_TRUE(r.contains(10));  // Inclusive start
    EXPECT_TRUE(r.contains(15));
    EXPECT_TRUE(r.contains(19));
    EXPECT_FALSE(r.contains(20)); // Exclusive stop
}

TEST(RangeTest, OperatorBracket_ReturnsOffsetValue) {
    Range r(100, 200);
    EXPECT_EQ(r[0], 100u);
    EXPECT_EQ(r[10], 110u);
    EXPECT_EQ(r[99], 199u);
}



// --------------------------
// Iterator & STL Integration
// --------------------------

TEST(RangeTest, Iterator_ForwardIteration) {
    Range r(10, 13);
    auto it = r.begin();

    EXPECT_EQ(*it, 10u);
    ++it;
    EXPECT_EQ(*it, 11u);
    it++;
    EXPECT_EQ(*it, 12u);
    ++it;
    EXPECT_EQ(it, r.end());
}

TEST(RangeTest, Iterator_RandomAccessArithmetic) {
    Range r(0, 100);
    auto it = r.begin();

    // +n
    EXPECT_EQ(*(it + 5), 5u);
    EXPECT_EQ(*(5 + it), 5u);

    // +=
    it += 10;
    EXPECT_EQ(*it, 10u);

    // -=
    it -= 5;
    EXPECT_EQ(*it, 5u);

    // -n
    EXPECT_EQ(*(it - 2), 3u);

    // difference
    auto it2 = r.begin() + 20;
    EXPECT_EQ(it2 - it, 15); // 20 - 5
}

TEST(RangeTest, Iterator_Comparison) {
    Range r(0, 10);
    auto it1 = r.begin();
    auto it2 = r.begin() + 5;

    EXPECT_LT(it1, it2);
    EXPECT_LE(it1, it2);
    EXPECT_GT(it2, it1);
    EXPECT_NE(it1, it2);

    it1 += 5;
    EXPECT_EQ(it1, it2);
}

TEST(RangeTest, WorksWithStdAlgorithms) {
    Range r(1, 6); // 1, 2, 3, 4, 5

    // std::find
    auto it = std::find(r.begin(), r.end(), 3);
    EXPECT_NE(it, r.end());
    EXPECT_EQ(*it, 3u);

    // std::accumulate
    size_t sum = std::accumulate(r.begin(), r.end(), 0ull);
    EXPECT_EQ(sum, 15u); // 1+2+3+4+5

    // std::ranges::copy
    std::vector<size_t> vec;
    std::ranges::copy(r, std::back_inserter(vec));

    ASSERT_EQ(vec.size(), 5u);
    EXPECT_EQ(vec[0], 1u);
    EXPECT_EQ(vec[4], 5u);
}

TEST(RangeTest, ConceptChecks) {
    static_assert(std::ranges::range<Range>);
    static_assert(std::ranges::sized_range<Range>);
    static_assert(std::ranges::random_access_range<Range>);
    static_assert(std::ranges::common_range<Range>); // begin/end same type
}











