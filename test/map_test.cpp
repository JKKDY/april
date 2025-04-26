#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <stdexcept>

#include <april/utils/map.hpp>



struct Dummy {
    int id;
};

using namespace april::utils::impl;

TEST(DenseMapTest, BuildThrowsOnSizeMismatch) {
    DensePairMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{1}));

    EXPECT_THROW(map.build(keys, std::move(values)), std::invalid_argument);
}

TEST(DenseMapTest, BuildThrowsOnDuplicateKeys) {
    DensePairMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {2,1} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{1}));
    values.emplace_back(std::make_unique<Dummy>(Dummy{2}));

    EXPECT_THROW(map.build(keys, std::move(values)), std::invalid_argument);
}

TEST(DenseMapTest, QueryElementPresent) {
    DensePairMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{100}));
    values.emplace_back(std::make_unique<Dummy>(Dummy{200}));

    map.build(keys, std::move(values));

    // Both orders should work
    Dummy* result1 = map.get(1,2);
    ASSERT_NE(result1, nullptr);
    EXPECT_EQ(result1->id, 100);

    Dummy* result2 = map.get(2,1);
    ASSERT_NE(result2, nullptr);
    EXPECT_EQ(result2->id, 100);

    Dummy* result3 = map.get(3,4);
    ASSERT_NE(result3, nullptr);
    EXPECT_EQ(result3->id, 200);

    Dummy* result4 = map.get(4,3);
    ASSERT_NE(result4, nullptr);
    EXPECT_EQ(result4->id, 200);
}

TEST(DenseMapTest, QueryElementAbsent) {
   DensePairMap<Dummy> map;
   std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
   std::vector<std::unique_ptr<Dummy>> values;
   values.emplace_back(std::make_unique<Dummy>(Dummy{100}));
   values.emplace_back(std::make_unique<Dummy>(Dummy{200}));

   map.build(keys, std::move(values));

   // Query a non-existing pair
   EXPECT_EQ(map.get(0,0), nullptr);
   EXPECT_EQ(map.get(2,3), nullptr);
#ifndef NDEBUG
   EXPECT_DEATH(map.get(5,0), ".*out of range.*");
#endif
}



TEST(UnorderedMapTest, BuildThrowsOnSizeMismatch) {
    UnorderedMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{1}));

    EXPECT_THROW(map.build(keys, std::move(values)), std::invalid_argument);
}

TEST(UnorderedMapTest, BuildThrowsOnDuplicateKeys) {
    UnorderedMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {2,1} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{1}));
    values.emplace_back(std::make_unique<Dummy>(Dummy{2}));

    EXPECT_THROW(map.build(keys, std::move(values)), std::invalid_argument);
}

TEST(UnorderedMapTest, QueryElementPresent) {
    UnorderedMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{100}));
    values.emplace_back(std::make_unique<Dummy>(Dummy{200}));

    map.build(keys, std::move(values));

    // Both orders should work
    Dummy* result1 = map.get(1,2);
    ASSERT_NE(result1, nullptr);
    EXPECT_EQ(result1->id, 100);

    Dummy* result2 = map.get(2,1);
    ASSERT_NE(result2, nullptr);
    EXPECT_EQ(result2->id, 100);

    Dummy* result3 = map.get(3,4);
    ASSERT_NE(result3, nullptr);
    EXPECT_EQ(result3->id, 200);

    Dummy* result4 = map.get(4,3);
    ASSERT_NE(result4, nullptr);
    EXPECT_EQ(result4->id, 200);
}

TEST(UnorderedMapTest, QueryElementAbsent) {
    UnorderedMap<Dummy> map;
    std::vector<std::pair<size_t, size_t>> keys = { {1,2}, {3,4} };
    std::vector<std::unique_ptr<Dummy>> values;
    values.emplace_back(std::make_unique<Dummy>(Dummy{100}));
    values.emplace_back(std::make_unique<Dummy>(Dummy{200}));

    map.build(keys, std::move(values));

    // Query a non-existing pair
    EXPECT_EQ(map.get(0,0), nullptr);
    EXPECT_EQ(map.get(2,3), nullptr);
    EXPECT_EQ(map.get(5,0), nullptr);
}