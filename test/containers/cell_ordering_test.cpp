#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_set>


#include "april/common.hpp"
#include "april/containers/cell_orderings.hpp"

using namespace april;

// ----------------------------
// HELPER FUNCTIONS FOR TESTING
// ----------------------------

// Reconstructs (x,y,z) from a flat index for a specific grid size
uint3 Unflatten(uint32_t flat_index, const uint3& dims) {
    uint3 p;
    p.x = flat_index % dims.x;
    p.y = (flat_index / dims.x) % dims.y;
    p.z = flat_index / (dims.x * dims.y);
    return p;
}

// Manhattan distance between two points
uint32_t Dist(const uint3& a, const uint3& b) {
    // Avoid underflow with std::abs on signed diff or simple logic
    auto diff = [](uint32_t u, uint32_t v) { return (u > v) ? u - v : v - u; };
    return diff(a.x, b.x) + diff(a.y, b.y) + diff(a.z, b.z);
}

// --------------------
// HILBERT TEST FIXTURE
// --------------------

class HilbertTest : public testing::Test {
protected:
    // Helper to run full verification on a grid
    void VerifyGrid(uint3 dims) {

        const size_t N = dims.x * dims.y * dims.z;
        const std::vector<uint32_t> ranking = container::hilbert_order(dims);

        // 1. check Result Size
        ASSERT_EQ(ranking.size(), N) << "Result vector size mismatch";

        // 2. Verify Bijectivity (Uniqueness)
        // The output 'ranking' maps FlatIndex -> CurveRank. 
        // We need to invert it to walk the curve: CurveRank -> FlatIndex.
        std::vector<uint32_t> curve_path(N);
        std::vector<bool> seen(N, false);

        for (size_t flat_idx = 0; flat_idx < N; ++flat_idx) {
            uint32_t rank = ranking[flat_idx];
            
            ASSERT_LT(rank, N) << "Rank out of bounds";
            ASSERT_FALSE(seen[rank]) << "Duplicate rank found: " << rank;
            
            seen[rank] = true;
            curve_path[rank] = flat_idx;
        }

        // 3. Verify Monotonicity of Hilbert Keys (Correctness of Sort)
        // Calculate the raw Hilbert Key for every point in the sorted path
        // and ensure strictly increasing order.
        const uint32_t max_dim = std::max({dims.x, dims.y, dims.z});
        const int bits = std::bit_width(max_dim > 0 ? max_dim - 1 : 0);

        uint64_t prev_key = 0;
        bool first = true;

        for (const uint32_t flat_idx : curve_path) {
            const uint3 p = Unflatten(flat_idx, dims);
            uint64_t current_key = container::internal::hilbert_encode_3d(p.x, p.y, p.z, bits);

            if (!first) {
                // Key must be >= previous key. 
                // Since coordinates are unique, keys MUST be strictly greater 
                // unless we have a hash collision (unlikely with 64 bits for reasonable sizes)
                EXPECT_GT(current_key, prev_key) 
                    << "Sorting failed! Hilbert Keys not monotonic at Rank " 
                    << ranking[flat_idx];
            }
            prev_key = current_key;
            first = false;
        }

        // 3. Verify Locality (Adjacency)
        // For a perfect Power-of-Two cube, distance is ALWAYS 1.
        // For arbitrary shapes, the curve might jump over "empty" space in the bounding box.
        // We track the maximum jump.
        uint32_t max_jump = 0;

        for (size_t i = 0; i < N - 1; ++i) {
            uint3 p1 = Unflatten(curve_path[i], dims);
            uint3 p2 = Unflatten(curve_path[i+1], dims);
            uint32_t d = Dist(p1, p2);

            if (d > max_jump) max_jump = d;
        }

        // Output stats for debugging
        // printf("Grid %u %u %u: Max Jump = %u, Avg Jump = %.2f\n", dims.x, dims.y, dims.z, max_jump, total_dist / (N-1));

        // CRITICAL CHECK:
        // If dimensions are Power of Two, Max Jump MUST be 1.
        if (N > 1) {
            const bool is_pow2 = (std::popcount(dims.x) == 1) &&
                           (std::popcount(dims.y) == 1) &&
                           (std::popcount(dims.z) == 1) &&
                           (dims.x == dims.y) && (dims.y == dims.z);

            if (is_pow2) {
                // For power-of-two cubes, step MUST be exactly 1
                EXPECT_EQ(max_jump, 1) << "FAILED at " << dims.x << "x" << dims.y << "x" << dims.z;
            } else {
                // For arbitrary shapes, larger jumps are allowed
                EXPECT_LT(max_jump, 20);
            }
        }
    }
};

// ---------------------------------------------------------------------------
// TESTS
// ---------------------------------------------------------------------------

TEST_F(HilbertTest, TrivialCase_1x1x1) {
    VerifyGrid({1, 1, 1});
}

TEST_F(HilbertTest, SmallCube_2x2x2) {
    VerifyGrid({2, 2, 2});
}

TEST_F(HilbertTest, StandardCube_4x4x4) {
    // gold standard test: max jump must be 1
    VerifyGrid({4, 4, 4});
}

TEST_F(HilbertTest, StandardCube_8x8x8) {
    // gold standard test: max jump must be 1
    VerifyGrid({8, 8, 8});
}

TEST_F(HilbertTest, StandardCube_16x16x16) {
    // gold standard test: max jump must be 1
    VerifyGrid({16, 16, 16});
}

TEST_F(HilbertTest, StandardCube_32x32x32) {
    // gold standard test: max jump must be 1
    VerifyGrid({32, 32, 32});
}

TEST_F(HilbertTest, StandardCube_64) {
    // gold standard test: max jump must be 1
    VerifyGrid({64, 64, 64});
}

TEST_F(HilbertTest, Rectangular_PowerOfTwo_4x4x2) {
    // Non-cubic, but fits perfectly in octants.
    VerifyGrid({4, 4, 2});
}

TEST_F(HilbertTest, NonPowerOfTwo_3x3x3) {
    // This tests the "Virtual Bounding Box" logic.
    // 3x3x3 fits in 4x4x4.
    VerifyGrid({3, 3, 3});
}

TEST_F(HilbertTest, Flat_Plate_10x10x1) {
    // 2D case effectively
    VerifyGrid({10, 10, 1});
}

TEST_F(HilbertTest, Prime_Dimensions_13x7x5) {
    // Stress test for arbitrary odd sizes
    VerifyGrid({13, 7, 5});
}



// ---------------------------------------------------------------------------
// MORTON TEST FIXTURE
// ---------------------------------------------------------------------------

class MortonTest : public testing::Test {
protected:
    // Helper to verify the bit-splitting logic specifically
    // 0b111 (7) should become 0b001001001 (73)
    // 0b11 (3)  should become 0b001001 (9)
    void VerifyBitSplit() {
        EXPECT_EQ(container::internal::split_by_3(0), 0);
        EXPECT_EQ(container::internal::split_by_3(1), 1);
        EXPECT_EQ(container::internal::split_by_3(2), 8);  // 10 -> 001000
        EXPECT_EQ(container::internal::split_by_3(3), 9);  // 11 -> 001001
        EXPECT_EQ(container::internal::split_by_3(7), 73); // 111 -> 001001001

        // Check max valid input (21 bits) preserves high bit
        uint32_t max_input = (1 << 21) - 1;
        uint64_t max_output = container::internal::split_by_3(max_input);
        EXPECT_GT(max_output, 0);
    }

    void VerifyGrid(uint3 dims) {
        const size_t N = dims.x * dims.y * dims.z;
        const std::vector<uint32_t> ranking = container::morton_order(dims);

        ASSERT_EQ(ranking.size(), N) << "Result vector size mismatch";

        // 1. Bijectivity Check
        std::vector<uint32_t> curve_path(N);
        std::vector<bool> seen(N, false);

        for (size_t flat_idx = 0; flat_idx < N; ++flat_idx) {
            uint32_t rank = ranking[flat_idx];
            ASSERT_LT(rank, N);
            ASSERT_FALSE(seen[rank]) << "Duplicate rank found: " << rank;

            seen[rank] = true;
            curve_path[rank] = flat_idx; // Store inverted mapping
        }

        // 2. Monotonicity Check (The Core Morton Property)
        // If we follow the path 0 -> 1 -> 2 ... -> N, the computed Morton codes
        // MUST be strictly increasing.
        uint64_t prev_key = 0;
        bool first = true;

        for (uint32_t flat_idx : curve_path) {
            const uint3 p = Unflatten(flat_idx, dims);
            uint64_t current_key = container::internal::morton_3d_64(p.x, p.y, p.z);

            if (!first) {
                // If this fails, the sorting logic inside morton_order is wrong
                EXPECT_GT(current_key, prev_key)
                    << "Morton keys not monotonic at Rank " << ranking[flat_idx]
                    << " for point " << p.x << "," << p.y << "," << p.z;
            }
            prev_key = current_key;
            first = false;
        }
    }
};

// -------------------
// MORTON TEST FIXTURE
// -------------------

TEST_F(MortonTest, BitManipulationCorrectness) {
    VerifyBitSplit();
}

TEST_F(MortonTest, Trivial_1x1x1) {
    VerifyGrid({1, 1, 1});
}

TEST_F(MortonTest, SimpleCube_2x2x2) {
    VerifyGrid({2, 2, 2});

    // Manual verification for 2x2x2 to ensure it's actually Z-order
    // Flat indices:
    // z=0: 0(0,0), 1(1,0), 2(0,1), 3(1,1)
    // z=1: 4(0,0), 5(1,0), 6(0,1), 7(1,1) (offset by 4)

    // Expected Z-Curve path of flat indices: 0, 1, 2, 3, 4, 5, 6, 7
    // Because Z-order on 2x2x2 is identical to linear order if x is fastest axis
    // Let's verify our specific mapping.

    // If our loops are:
    // for z { for y { for x ... } }
    // Then Linear Index IS Morton Index for 2x2x2.
    auto result = container::morton_order({2,2,2});
    for(size_t i=0; i<8; ++i) {
        EXPECT_EQ(result[i], i) << "For 2x2x2, Morton order should match Linear order";
    }
}

TEST_F(MortonTest, StandardCube_4x4x4) {
    VerifyGrid({4, 4, 4});
}

TEST_F(MortonTest, StandardCube_8x8x8) {
    VerifyGrid({8, 8, 8});
}

TEST_F(MortonTest, StandardCube_32x32x32) {
    VerifyGrid({32, 32, 32});
}

TEST_F(MortonTest, Rectangular_4x4x2) {
    VerifyGrid({4, 4, 2});
}

TEST_F(MortonTest, NonPowerOfTwo_3x3x3) {
    // Morton codes work fine for non-powers of two,
    // they just skip the codes for the missing coordinates.
    VerifyGrid({3, 3, 3});
}

TEST_F(MortonTest, Flat_Plate_10x10x1) {
    VerifyGrid({10, 10, 1});
}

TEST_F(MortonTest, Prime_Dimensions_13x7x5) {
    VerifyGrid({13, 7, 5});
}