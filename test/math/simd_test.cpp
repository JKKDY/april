#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "april/simd/backend_std_simd.hpp"
#include "april/simd/backend_xsimd.hpp"
#include "april/simd/concepts.hpp"

using BackendTypes = testing::Types<
    april::simd::internal::xsimd::Wide<double>,
    april::simd::internal::std_simd::Wide<double>,
    april::simd::internal::xsimd::Wide<float>,
    april::simd::internal::std_simd::Wide<float>
>;

AP_SIMD_IMPORT_WIDE_MATH(xsimd)
AP_SIMD_IMPORT_WIDE_MATH(std_simd)

template <typename T>
class SimdWideTest : public testing::Test {
public:
    using WideType = T;
    using Scalar = T::value_type;
    static constexpr size_t Size = T::size();
};

TYPED_TEST_SUITE(SimdWideTest, BackendTypes);

// Concept Verification
TYPED_TEST(SimdWideTest, SatisfiesSimdConcept) {
    static_assert(april::simd::IsSimdType<typename TestFixture::WideType>,
        "Type does not satisfy IsSimdType concept");
}

// Load / Store / Broadcast
TYPED_TEST(SimdWideTest, LoadStoreBroadcast) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    // Broadcast
    Wide w_scalar(42.0);
    std::vector<Scalar> buffer(N);
    w_scalar.store(buffer.data());

    for (auto v : buffer) {
        EXPECT_DOUBLE_EQ(v, 42.0);
    }

    // Load
    std::iota(buffer.begin(), buffer.end(), 0.0); // 0, 1, 2...
    Wide w_load = Wide::load(buffer.data());

    std::vector<Scalar> out(N);
    w_load.store(out.data());

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(out[i], static_cast<Scalar>(i));
    }
}

// Check Arithmetic Operations
TYPED_TEST(SimdWideTest, Arithmetic) {
    using Wide = typename TestFixture::WideType;
    using Scalar = typename TestFixture::Scalar;
    size_t N = TestFixture::Size;

    Wide a(10.0);
    Wide b(2.0);

    // Basic Ops
    Wide sum = a + b;
    Wide diff = a - b;
    Wide prod = a * b;
    Wide quot = a / b;

    std::vector<Scalar> res_sum(N), res_diff(N), res_prod(N), res_quot(N);
    sum.store(res_sum.data());
    diff.store(res_diff.data());
    prod.store(res_prod.data());
    quot.store(res_quot.data());

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(res_sum[i], 12.0);
        EXPECT_DOUBLE_EQ(res_diff[i], 8.0);
        EXPECT_DOUBLE_EQ(res_prod[i], 20.0);
        EXPECT_DOUBLE_EQ(res_quot[i], 5.0);
    }

    // Compound Assignment
    a += b; // a is now 12
    std::vector<Scalar> res_compound(N);
    a.store(res_compound.data());
    for (auto v : res_compound) EXPECT_DOUBLE_EQ(v, 12.0);
}

// Check free Math Functions (ADL Check)
TYPED_TEST(SimdWideTest, MathFunctions) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    Wide val(16.0);

    // Using ADL (Argument Dependent Lookup) - No explicit namespace needed
    // The "Hidden Friend" functions inside the struct enable this.
    Wide res_sqrt = sqrt(val);
    Wide res_rsqrt = rsqrt(val);

    std::vector<Scalar> s_out(N), r_out(N);
    res_sqrt.store(s_out.data());
    res_rsqrt.store(r_out.data());

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(s_out[i], 4.0);
        EXPECT_NEAR(r_out[i], 0.25, 1e-4);
    }
}

// Rotation
TYPED_TEST(SimdWideTest, Rotation) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    // Setup: [0, 1, 2, 3 ...]
    std::vector<Scalar> data(N);
    std::iota(data.begin(), data.end(), 0.0);

    Wide w = Wide::load(data.data());

    // Rotate Left by 1 -> [1, 2, 3, 0]
    Wide w_rot = w.rotate_left();

    std::vector<Scalar> out(N);
    w_rot.store(out.data());

    for (size_t i = 0; i < N; ++i) {
        Scalar expected = data[(i + 1) % N];
        EXPECT_DOUBLE_EQ(out[i], expected) << "Mismatch at index " << i;
    }
}

// Check Gather (Indirect Load)
TYPED_TEST(SimdWideTest, Gather) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    // Source memory: [100, 101, 102, ... 120]
    std::vector<Scalar> memory(100);
    std::iota(memory.begin(), memory.end(), 100.0);

    // Create pointer array: Pick indices 0, 2, 4, 6...
    std::vector<const Scalar*> ptrs(N);
    for(size_t i=0; i<N; ++i) {
        ptrs[i] = &memory[i * 2];
    }

    // Perform Gather
    // Note: We use the pointer-array gather signature
    Wide gathered = Wide::gather(ptrs.data());

    std::vector<Scalar> out(N);
    gathered.store(out.data());

    for (size_t i = 0; i < N; ++i) {
        double expected = 100.0 + (i * 2.0);
        EXPECT_DOUBLE_EQ(out[i], expected);
    }
}