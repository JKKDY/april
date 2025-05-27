#include <gtest/gtest.h>

#include "april/utils/set.hpp"

// For death tests
#ifdef NDEBUG
#  define MAYBE_DEATH_TEST(test_suite_name, test_name) DISABLED_##test_name
#else
#  define MAYBE_DEATH_TEST(test_suite_name, test_name) test_name
#endif

using namespace april::utils;
using UInt = uint32_t;

TEST(IndexSet, ConstructEmpty) {
    IndexSet<UInt> s{9};         // allows ids 0..9
    EXPECT_EQ(s.size(), 0u);
    // nothing contained
    EXPECT_FALSE(s.contains(0));
    EXPECT_FALSE(s.contains(9));
    EXPECT_FALSE(s.contains(10));  // out of range
}

TEST(IndexSet, SingleInsertContainsAndSize) {
    IndexSet<UInt> s{5};
    s.insert(3);
    EXPECT_TRUE(s.contains(3));
    EXPECT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0], 3u);
    // other ids still absent
    EXPECT_FALSE(s.contains(2));
}

TEST(IndexSet, MultipleInsertsContainsAll) {
    IndexSet<UInt> s{100};
    s.insert(10);
    s.insert(42);
    s.insert(7);
    EXPECT_EQ(s.size(), 3u);

    // collect
    std::vector seen(s.begin(), s.end());
    std::sort(seen.begin(), seen.end());

    const std::vector<UInt> expected{7,10,42};
    EXPECT_EQ(seen, expected);
}


TEST(IndexSet, EraseRemovesAndSwapsBack) {
    IndexSet<UInt> s{10};
    s.insert(2);
    s.insert(5);
    s.insert(8);
    ASSERT_EQ(s.size(), 3u);

    // Erase middle
    s.erase(5);
    EXPECT_FALSE(s.contains(5));
    EXPECT_EQ(s.size(), 2u);

    // The element 8 should have been swapped into position of 5
    std::vector remaining(s.begin(), s.end());
    std::sort(remaining.begin(), remaining.end());
    EXPECT_EQ(remaining, std::vector<UInt>({2,8}));
}

TEST(IndexSet, ReinsertAfterErase) {
    IndexSet<UInt> s{3};
    s.insert(1);
    s.erase(1);
    EXPECT_FALSE(s.contains(1));
    EXPECT_EQ(s.size(), 0u);

    // Now we can re-insert
    s.insert(1);
    EXPECT_TRUE(s.contains(1));
    EXPECT_EQ(s.size(), 1u);
    EXPECT_EQ(s[0], 1u);
}

// Death tests for invalid operations
TEST(MAYBE_DEATH_TEST(IndexSet, InsertDuplicate), IndexSetError_InsertDuplicate) {
    IndexSet<UInt> s{2};
    s.insert(0);
    EXPECT_DEATH(s.insert(0), "");  // duplicate
}

TEST(MAYBE_DEATH_TEST(IndexSet, InsertOutOfRange), IndexSetError_InsertOutOfRange) {
    IndexSet<UInt> s{2};
    EXPECT_DEATH(s.insert(3), "");  // 3 >= N=3
}

TEST(MAYBE_DEATH_TEST(IndexSet, EraseNonexistent), IndexSetError_EraseNonexistent) {
    IndexSet<UInt> s{5};
    EXPECT_DEATH(s.erase(1), "");  // never inserted
}

TEST(MAYBE_DEATH_TEST(IndexSet, EraseOutOfRange), IndexSetError_EraseOutOfRange) {
    IndexSet<UInt> s{5};
    EXPECT_DEATH(s.erase(7), "");  // 7 >= N
}

TEST(IndexSet, OutOfRangeReturnsFalse) {
    IndexSet<UInt> s{4};
    // contains should simply return false for out-of-range
    EXPECT_FALSE(s.contains(10u));
}

// Test for iterator begin/end on empty
TEST(IndexSet, EmptyBeginEqualsEnd) {
    IndexSet<UInt> s{0};
    EXPECT_EQ(s.begin(), s.end());
}

// Test large sequence
TEST(IndexSet, ManyInsertsAndErases) {
    const UInt MAX_ID = 1000;
    IndexSet<UInt> s{MAX_ID};
    // insert all even IDs
    for (UInt i = 0; i <= MAX_ID; i += 2) {
        s.insert(i);
    }
    EXPECT_EQ(s.size(), (MAX_ID/2) + 1);
    for (UInt i = 0; i <= MAX_ID; ++i) {
        EXPECT_EQ(s.contains(i), (i % 2 == 0));
    }

    // erase all evens
    for (UInt i = 0; i <= MAX_ID; i += 2) {
        s.erase(i);
    }
    EXPECT_EQ(s.size(), 0u);
    for (UInt i = 0; i <= MAX_ID; ++i) {
        EXPECT_FALSE(s.contains(i));
    }
}
