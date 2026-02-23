#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

#include "april/simd/backend_std_simd.hpp"
#include "april/simd/backend_xsimd.hpp"
#include "april/simd/simd_traits.hpp"

using BackendTypes = testing::Types<
    april::simd::internal::xsimd::Packed<double>,
    april::simd::internal::std_simd::Packed<double>,
    april::simd::internal::xsimd::Packed<float>,
    april::simd::internal::std_simd::Packed<float>
>;


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
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
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

// check mixed arithmetic (wide + scalar)
TYPED_TEST(SimdWideTest, ArithmeticScalarLeftf) {
    using Wide = TestFixture::WideType;
    size_t N = TestFixture::Size;

    float a = 10.0;
    Wide b = 2.0;

    // Basic Ops
    Wide sum = a + b;
    Wide diff = a - b;
    Wide prod = a * b;
    Wide quot = a / b;

    auto res_sum =  sum.to_array();
    auto res_diff = diff.to_array();
    auto res_prod = prod.to_array();
    auto res_quot = quot.to_array();

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(res_sum [i], 12.0);
        EXPECT_DOUBLE_EQ(res_diff[i], 8.0);
        EXPECT_DOUBLE_EQ(res_prod[i], 20.0);
        EXPECT_DOUBLE_EQ(res_quot[i], 5.0);
    }
}

TYPED_TEST(SimdWideTest, ArithmeticScalarLeftd) {
    using Wide = TestFixture::WideType;
    size_t N = TestFixture::Size;

    double a = 10.0;
    Wide b = 2.0;

    // Basic Ops
    Wide sum = a + b;
    Wide diff = a - b;
    Wide prod = a * b;
    Wide quot = a / b;

    auto res_sum =  sum.to_array();
    auto res_diff = diff.to_array();
    auto res_prod = prod.to_array();
    auto res_quot = quot.to_array();

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(res_sum [i], 12.0);
        EXPECT_DOUBLE_EQ(res_diff[i], 8.0);
        EXPECT_DOUBLE_EQ(res_prod[i], 20.0);
        EXPECT_DOUBLE_EQ(res_quot[i], 5.0);
    }
}

TYPED_TEST(SimdWideTest, ArithmeticScalarRightF) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    Wide a = 10.0;
    float b = 2.0;

    // Basic Ops
    Wide sum = a + b;
    Wide diff = a - b;
    Wide prod = a * b;
    Wide quot = a / b;

    auto res_sum =  sum.to_array();
    auto res_diff = diff.to_array();
    auto res_prod = prod.to_array();
    auto res_quot = quot.to_array();

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

TYPED_TEST(SimdWideTest, ArithmeticScalarRightD) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    Wide a = 10.0;
    double b = 2.0;

    // Basic Ops
    Wide sum = a + b;
    Wide diff = a - b;
    Wide prod = a * b;
    Wide quot = a / b;

    auto res_sum =  sum.to_array();
    auto res_diff = diff.to_array();
    auto res_prod = prod.to_array();
    auto res_quot = quot.to_array();

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


// Check Scalar Assignment (Broadcasting)
TYPED_TEST(SimdWideTest, ScalarAssignment) {
    using Wide = TestFixture::WideType;
    const size_t N = TestFixture::Size;

    Wide w;
    w = 7.0; // Trigger the assignment operator we added

    auto out = w.to_array();
    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(out[i], 7.0);
    }
}

// Check Conditional Select and Masking
TYPED_TEST(SimdWideTest, SelectAndMasking) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    std::vector<Scalar> vals_a(N), vals_b(N);
    for(size_t i = 0; i < N; ++i) {
        vals_a[i] = static_cast<Scalar>(i);       // 0, 1, 2, 3...
        vals_b[i] = static_cast<Scalar>(N - i);   // N, N-1, N-2...
    }

    Wide a = Wide::load(vals_a.data());
    Wide b = Wide::load(vals_b.data());

    // Create a mask (e.g., true where a < b)
    auto mask = a < b;

    // Select: result[i] = mask[i] ? a[i] : b[i]
    // This effectively computes the element-wise minimum.
    Wide res = select(mask, a, b);
    auto out = res.to_array();

    for (size_t i = 0; i < N; ++i) {
        Scalar expected = (vals_a[i] < vals_b[i]) ? vals_a[i] : vals_b[i];
        EXPECT_DOUBLE_EQ(out[i], expected);
    }
}

// Check Horizontal Reductions (Math)
TYPED_TEST(SimdWideTest, MathReductions) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    std::vector<Scalar> vals(N);
    Scalar expected_sum = 0;
    Scalar expected_min = std::numeric_limits<Scalar>::max();
    Scalar expected_max = std::numeric_limits<Scalar>::lowest();

    for(size_t i = 0; i < N; ++i) {
        vals[i] = static_cast<Scalar>(i + 1); // 1, 2, 3...
        expected_sum += vals[i];
        expected_min = std::min(expected_min, vals[i]);
        expected_max = std::max(expected_max, vals[i]);
    }

    Wide w = Wide::load(vals.data());

    EXPECT_DOUBLE_EQ(w.reduce_add(), expected_sum);
    EXPECT_DOUBLE_EQ(w.reduce_min(), expected_min);
    EXPECT_DOUBLE_EQ(w.reduce_max(), expected_max);
}

// Check Logical Mask Reductions
TYPED_TEST(SimdWideTest, LogicalMaskReductions) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;

    Wide zeros = 0.0;
    Wide ones = 1.0;

    auto mask_all_true = (zeros == zeros);
    auto mask_all_false = (zeros == ones);

    // Test homogeneous masks
    EXPECT_TRUE(all(mask_all_true));
    EXPECT_TRUE(any(mask_all_true));
    EXPECT_FALSE(none(mask_all_true));

    EXPECT_FALSE(all(mask_all_false));
    EXPECT_FALSE(any(mask_all_false));
    EXPECT_TRUE(none(mask_all_false));

    // Test heterogeneous mask (if vector size > 1)
    if constexpr (TestFixture::Size > 1) {
        std::vector<Scalar> vals(TestFixture::Size, 0.0);
        vals[0] = 1.0; // Make only the first element match

        Wide mix_val = Wide::load(vals.data());
        auto mask_mixed = (mix_val == ones);

        EXPECT_FALSE(all(mask_mixed));
        EXPECT_TRUE(any(mask_mixed));
        EXPECT_FALSE(none(mask_mixed));
    }
}








