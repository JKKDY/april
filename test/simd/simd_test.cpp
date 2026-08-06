#include <gtest/gtest.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <array>
#include <concepts>
#include <cstdint>
#include <type_traits>
#include <utility>


#include "april/simd/packed.hpp"
#include "april/simd/packed_concept.hpp"
#include "april/simd/backends/backend_scalar.hpp"



using BackendTypes = testing::Types<
    april::simd::Packed<double>,
    april::simd::Packed<float>,
    april::simd::internal::scalar::Packed<double>,
    april::simd::internal::scalar::Packed<float>
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
    std::iota(buffer.begin(), buffer.end(), static_cast<Scalar>(0)); // 0, 1, 2...
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
    const size_t N = TestFixture::Size;

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
    Wide res_sqrt = april::sqrt(val);
    Wide res_rsqrt = april::rsqrt(val);

    std::vector<Scalar> s_out(N), r_out(N);
    res_sqrt.store(s_out.data());
    res_rsqrt.store(r_out.data());

    for (size_t i = 0; i < N; ++i) {
        EXPECT_DOUBLE_EQ(s_out[i], 4.0);
        EXPECT_NEAR(r_out[i], 0.25, 1e-3);
    }
}

// Rotation
TYPED_TEST(SimdWideTest, Rotation) {
    using Wide = TestFixture::WideType;
    using Scalar = TestFixture::Scalar;
    size_t N = TestFixture::Size;

    // Setup: [0, 1, 2, 3 ...]
    std::vector<Scalar> data(N);
    std::iota(data.begin(), data.end(), Scalar(0.0));

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
    std::iota(memory.begin(), memory.end(), Scalar(100.0));

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
    Wide res = april::select(mask, a, b);
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
    EXPECT_TRUE(april::all(mask_all_true));
    EXPECT_TRUE(april::any(mask_all_true));
    EXPECT_FALSE(april::none(mask_all_true));

    EXPECT_FALSE(april::all(mask_all_false));
    EXPECT_FALSE(april::any(mask_all_false));
    EXPECT_TRUE(april::none(mask_all_false));

    // Test heterogeneous mask (if vector size > 1)
    if constexpr (TestFixture::Size > 1) {
        std::vector<Scalar> vals(TestFixture::Size, 0.0);
        vals[0] = 1.0; // Make only the first element match

        Wide mix_val = Wide::load(vals.data());
        auto mask_mixed = (mix_val == ones);

        EXPECT_FALSE(april::all(mask_mixed));
        EXPECT_TRUE(april::any(mask_mixed));
        EXPECT_FALSE(april::none(mask_mixed));
    }
}



// ---------------------------------------------------------
// NARROW TYPES TEST SUITE (Upcast/Downcast Memory Interface)
// ---------------------------------------------------------

using NarrowBackendTypes = testing::Types<
    april::simd::Packed<uint64_t>,
    april::simd::Packed<uint32_t>
>;

template <typename T>
class SimdNarrowTest : public testing::Test {
public:
    using WideType = T;
    using WideScalar = T::value_type;
    static constexpr size_t Size = T::size();
};

TYPED_TEST_SUITE(SimdNarrowTest, NarrowBackendTypes);

// Test 1: Load 8-bit integers into wide registers
TYPED_TEST(SimdNarrowTest, LoadNarrowUint8) {
    using Wide = TestFixture::WideType;
    using WideScalar = TestFixture::WideScalar;
    size_t N = TestFixture::Size;

    // Allocate EXACTLY N bytes to fiercely test for ASAN over-reads
    std::vector<uint8_t> narrow_mem(N);
    for (size_t i = 0; i < N; ++i) {
        narrow_mem[i] = static_cast<uint8_t>(i + 10);
    }

    Wide w_load = Wide::load(narrow_mem.data());

    // Export to wide array to check upcasting
    auto out = w_load.to_array();

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(out[i], static_cast<WideScalar>(i + 10));
    }
}

// Test 2: Store wide registers into 8-bit memory (Truncation)
TYPED_TEST(SimdNarrowTest, StoreNarrowUint8) {
    using Wide = TestFixture::WideType;
    size_t N = TestFixture::Size;

    // Setup wide data: [20, 21, 22...]
    std::vector<typename TestFixture::WideScalar> wide_mem(N);
    std::iota(wide_mem.begin(), wide_mem.end(), 20);
    Wide w_val = Wide::load(wide_mem.data());

    // Allocate EXACTLY N bytes to test for ASAN over-writes
    std::vector<uint8_t> narrow_mem(N, 0);
    w_val.store(narrow_mem.data());

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(narrow_mem[i], static_cast<uint8_t>(i + 20));
    }
}

// Test 3: Load 16-bit integers into wide registers
TYPED_TEST(SimdNarrowTest, LoadNarrowUint16) {
    using Wide = TestFixture::WideType;
    using WideScalar = TestFixture::WideScalar;
    size_t N = TestFixture::Size;

    std::vector<uint16_t> narrow_mem(N);
    for (size_t i = 0; i < N; ++i) {
        narrow_mem[i] = static_cast<uint16_t>(i + 1000);
    }

    Wide w_load = Wide::load(narrow_mem.data());
    auto out = w_load.to_array();

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(out[i], static_cast<WideScalar>(i + 1000));
    }
}

// Test 4: Gather narrow integers from non-contiguous memory
TYPED_TEST(SimdNarrowTest, GatherNarrow) {
    using Wide = TestFixture::WideType;
    using WideScalar = TestFixture::WideScalar;
    size_t N = TestFixture::Size;

    // Source memory: [50, 51, 52...]
    std::vector<uint8_t> memory(100);
    std::iota(memory.begin(), memory.end(), static_cast<uint8_t>(50));

    // Pointer array: Pick indices 0, 2, 4, 6...
    std::vector<const uint8_t*> ptrs(N);
    for(size_t i = 0; i < N; ++i) {
        ptrs[i] = &memory[i * 2];
    }

    Wide gathered = Wide::gather(ptrs.data());
    auto out = gathered.to_array();

    for (size_t i = 0; i < N; ++i) {
        EXPECT_EQ(out[i], static_cast<WideScalar>(50 + (i * 2)));
    }
}



// ---------------------------------------------------------
// EXPLICIT WIDTH AND PACKED/MASK INTEROPERABILITY TESTS
// ---------------------------------------------------------

namespace {

    template<typename PackedT>
    using PackedMaskOf = decltype(
        std::declval<const PackedT&>() ==
        std::declval<const PackedT&>()
    );


    template<typename ScalarT, size_t WidthV>
    struct SimdWidthCase {
        using Scalar = ScalarT;
        using PackedT = april::simd::Packed<Scalar, WidthV>;

        static constexpr size_t Width = WidthV;
    };


    template<typename MaskT>
    MaskT MakeAlternatingMask(bool first_lane_active = true) {
        std::array<bool, MaskT::size()> lanes{};

        for (size_t i = 0; i < MaskT::size(); ++i) {
            lanes[i] = ((i % 2) == 0) == first_lane_active;
        }

        return MaskT::load_unaligned(lanes.data());
    }


    template<typename MaskT>
    MaskT MakeMaskPattern() {
        std::array<bool, MaskT::size()> lanes{};

        for (size_t i = 0; i < MaskT::size(); ++i) {
            lanes[i] = (i % 3) != 1;
        }

        return MaskT::load_unaligned(lanes.data());
    }


    template<typename MaskA, typename MaskB>
    void ExpectSameMask(const MaskA& lhs, const MaskB& rhs) {
        static_assert(MaskA::size() == MaskB::size());

        const auto lhs_lanes = lhs.to_array();
        const auto rhs_lanes = rhs.to_array();

        for (size_t i = 0; i < MaskA::size(); ++i) {
            EXPECT_EQ(lhs_lanes[i], rhs_lanes[i])
                << "Mask mismatch at lane " << i;
        }
    }


    template<typename PackedT>
    void ExpectPackedValues(
        const PackedT& packed,
        const std::array<typename PackedT::value_type, PackedT::size()>& expected
    ) {
        const auto actual = packed.to_array();

        for (size_t i = 0; i < PackedT::size(); ++i) {
            EXPECT_EQ(actual[i], expected[i])
                << "Packed mismatch at lane " << i;
        }
    }

}


// Widths two and four should be available on the baseline SIMD targets
// already required by these packed types.
using ExplicitWidthCases = testing::Types<
    SimdWidthCase<float, 4>,
    SimdWidthCase<double, 2>,
    SimdWidthCase<uint32_t, 4>,
    SimdWidthCase<uint64_t, 2>
>;


template<typename T>
class SimdExplicitWidthTest : public testing::Test {
public:
    using Case = T;
    using PackedT =  Case::PackedT;
    using Scalar = Case::Scalar;
    using MaskT = PackedMaskOf<PackedT>;

    static constexpr size_t Width = Case::Width;
};


TYPED_TEST_SUITE(SimdExplicitWidthTest, ExplicitWidthCases);


// ---------------------------------------------------------
// EXPLICIT WIDTH TYPE CONTRACT
// ---------------------------------------------------------

TYPED_TEST(SimdExplicitWidthTest, ReportsRequestedWidth) {
    using PackedT = typename TestFixture::PackedT;
    using MaskT = typename TestFixture::MaskT;

    static_assert(april::simd::IsSimdType<PackedT>);
    static_assert(april::simd::IsSimdMask<MaskT>);

    static_assert(PackedT::size() == TestFixture::Width);
    static_assert(MaskT::size() == TestFixture::Width);

    EXPECT_EQ(PackedT::size(), TestFixture::Width);
    EXPECT_EQ(MaskT::size(), TestFixture::Width);
}


TYPED_TEST(SimdExplicitWidthTest, ComparisonsReturnMatchingWidthMask) {
    using PackedT = typename TestFixture::PackedT;
    using MaskT = typename TestFixture::MaskT;

    static_assert(std::same_as<
        decltype(std::declval<PackedT>() == std::declval<PackedT>()),
        MaskT
    >);

    static_assert(std::same_as<
        decltype(std::declval<PackedT>() < std::declval<PackedT>()),
        MaskT
    >);
}


// ---------------------------------------------------------
// EXPLICIT WIDTH LOAD/STORE
// ---------------------------------------------------------

TYPED_TEST(SimdExplicitWidthTest, LoadStorePreservesEveryLane) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;

    std::array<Scalar, TestFixture::Width> input{};

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        input[i] = static_cast<Scalar>(i + 10);
    }

    PackedT packed = PackedT::load_unaligned(input.data());

    std::array<Scalar, TestFixture::Width> output{};
    packed.store_unaligned(output.data());

    EXPECT_EQ(input, output);
}


TYPED_TEST(SimdExplicitWidthTest, BroadcastUsesRequestedWidth) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;

    PackedT packed(Scalar{7});
    const auto lanes = packed.to_array();

    static_assert(std::tuple_size_v<decltype(lanes)> == TestFixture::Width);

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        EXPECT_EQ(lanes[i], Scalar{7});
    }
}


// ---------------------------------------------------------
// EXPLICIT WIDTH ARITHMETIC
// ---------------------------------------------------------

TYPED_TEST(SimdExplicitWidthTest, ArithmeticOperatesOnRequestedWidth) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;

    PackedT lhs(Scalar{10});
    PackedT rhs(Scalar{2});

    const auto sum = (lhs + rhs).to_array();
    const auto difference = (lhs - rhs).to_array();
    const auto product = (lhs * rhs).to_array();
    const auto quotient = (lhs / rhs).to_array();

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        EXPECT_EQ(sum[i], Scalar{12});
        EXPECT_EQ(difference[i], Scalar{8});
        EXPECT_EQ(product[i], Scalar{20});
        EXPECT_EQ(quotient[i], Scalar{5});
    }
}


TYPED_TEST(SimdExplicitWidthTest, CompoundArithmeticOperatesOnRequestedWidth) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;

    PackedT value(Scalar{10});
    PackedT rhs(Scalar{2});

    value += rhs;
    value *= rhs;
    value -= rhs;
    value /= rhs;

    const auto lanes = value.to_array();

    // ((10 + 2) * 2 - 2) / 2 = 11
    for (size_t i = 0; i < TestFixture::Width; ++i) {
        EXPECT_EQ(lanes[i], Scalar{11});
    }
}


// ---------------------------------------------------------
// EXPLICIT WIDTH COMPARISON AND SELECT
// ---------------------------------------------------------

TYPED_TEST(SimdExplicitWidthTest, ComparisonProducesCorrectMaskPattern) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;

    std::array<Scalar, TestFixture::Width> lhs_values{};
    std::array<Scalar, TestFixture::Width> rhs_values{};
    std::array<bool, TestFixture::Width> expected{};

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        lhs_values[i] = static_cast<Scalar>(i);
        rhs_values[i] = static_cast<Scalar>(
            (i % 2) == 0 ? i + 1 : i
        );

        expected[i] = lhs_values[i] < rhs_values[i];
    }

    PackedT lhs = PackedT::load_unaligned(lhs_values.data());
    PackedT rhs = PackedT::load_unaligned(rhs_values.data());

    const auto actual = (lhs < rhs).to_array();

    EXPECT_EQ(actual, expected);
}


TYPED_TEST(SimdExplicitWidthTest, SelectUsesEveryMaskLane) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Width> mask_values{};
    std::array<Scalar, TestFixture::Width> true_values{};
    std::array<Scalar, TestFixture::Width> false_values{};
    std::array<Scalar, TestFixture::Width> expected{};

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        mask_values[i] = (i % 2) == 0;
        true_values[i] = static_cast<Scalar>(100 + i);
        false_values[i] = static_cast<Scalar>(200 + i);

        expected[i] = mask_values[i]
            ? true_values[i]
            : false_values[i];
    }

    MaskT mask = MaskT::load_unaligned(mask_values.data());
    PackedT true_value = PackedT::load_unaligned(true_values.data());
    PackedT false_value = PackedT::load_unaligned(false_values.data());

    PackedT result = select(mask, true_value, false_value);

    ExpectPackedValues(result, expected);
}


// ---------------------------------------------------------
// SAME-WIDTH, SAME-SIZE MASK CONVERSION
// ---------------------------------------------------------

TEST(SimdWidthInteroperabilityTest, Float4MaskConvertsToUint32Mask) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    static_assert(FloatPacked::size() == UIntPacked::size());
    static_assert(FloatMask::size() == UIntMask::size());
    static_assert(sizeof(float) == sizeof(uint32_t));
    static_assert(std::is_convertible_v<FloatMask, UIntMask>);

    FloatMask source = MakeMaskPattern<FloatMask>();
    UIntMask converted = static_cast<UIntMask>(source);

    ExpectSameMask(source, converted);
}


TEST(SimdWidthInteroperabilityTest, Uint32MaskConvertsToFloat4Mask) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    static_assert(std::is_convertible_v<UIntMask, FloatMask>);

    UIntMask source = MakeMaskPattern<UIntMask>();
    FloatMask converted = static_cast<FloatMask>(source);

    ExpectSameMask(source, converted);
}


TEST(SimdWidthInteroperabilityTest, Double2MaskConvertsToUint64Mask) {
    using DoublePacked = april::simd::Packed<double, 2>;
    using UIntPacked = april::simd::Packed<uint64_t, 2>;

    using DoubleMask = PackedMaskOf<DoublePacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    static_assert(DoublePacked::size() == UIntPacked::size());
    static_assert(DoubleMask::size() == UIntMask::size());
    static_assert(sizeof(double) == sizeof(uint64_t));
    static_assert(std::is_convertible_v<DoubleMask, UIntMask>);

    DoubleMask source = MakeMaskPattern<DoubleMask>();
    UIntMask converted = static_cast<UIntMask>(source);

    ExpectSameMask(source, converted);
}


TEST(SimdWidthInteroperabilityTest, Uint64MaskConvertsToDouble2Mask) {
    using DoublePacked = april::simd::Packed<double, 2>;
    using UIntPacked = april::simd::Packed<uint64_t, 2>;

    using DoubleMask = PackedMaskOf<DoublePacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    static_assert(std::is_convertible_v<UIntMask, DoubleMask>);

    UIntMask source = MakeMaskPattern<UIntMask>();
    DoubleMask converted = static_cast<DoubleMask>(source);

    ExpectSameMask(source, converted);
}


// ---------------------------------------------------------
// CONVERSION ROUND TRIPS
// ---------------------------------------------------------

TEST(SimdWidthInteroperabilityTest, FloatUint32MaskRoundTripPreservesLanes) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    FloatMask original = MakeAlternatingMask<FloatMask>();
    UIntMask integer_mask = static_cast<UIntMask>(original);
    FloatMask round_trip = static_cast<FloatMask>(integer_mask);

    ExpectSameMask(original, integer_mask);
    ExpectSameMask(original, round_trip);
}


TEST(SimdWidthInteroperabilityTest, DoubleUint64MaskRoundTripPreservesLanes) {
    using DoublePacked = april::simd::Packed<double, 2>;
    using UIntPacked = april::simd::Packed<uint64_t, 2>;

    using DoubleMask = PackedMaskOf<DoublePacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    DoubleMask original = MakeAlternatingMask<DoubleMask>(false);
    UIntMask integer_mask = static_cast<UIntMask>(original);
    DoubleMask round_trip = static_cast<DoubleMask>(integer_mask);

    ExpectSameMask(original, integer_mask);
    ExpectSameMask(original, round_trip);
}


// ---------------------------------------------------------
// MIXED-TYPE MASK OPERATORS
// ---------------------------------------------------------

TEST(SimdWidthInteroperabilityTest, FloatAndUint32MasksCanBeCombined) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    std::array<bool, FloatMask::size()> lhs_values{};
    std::array<bool, UIntMask::size()> rhs_values{};

    for (size_t i = 0; i < FloatMask::size(); ++i) {
        lhs_values[i] = (i % 2) == 0;
        rhs_values[i] = (i % 3) == 0;
    }

    FloatMask lhs = FloatMask::load_unaligned(lhs_values.data());
    UIntMask rhs = UIntMask::load_unaligned(rhs_values.data());

    auto result = lhs & rhs;

    static_assert(std::same_as<decltype(result), FloatMask>);

    const auto actual = result.to_array();

    for (size_t i = 0; i < FloatMask::size(); ++i) {
        EXPECT_EQ(actual[i], lhs_values[i] && rhs_values[i]);
    }
}


TEST(SimdWidthInteroperabilityTest, Uint32OrFloatMasksCanBeCombined) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    std::array<bool, FloatMask::size()> float_values{};
    std::array<bool, UIntMask::size()> uint_values{};

    for (size_t i = 0; i < FloatMask::size(); ++i) {
        float_values[i] = (i % 2) == 0;
        uint_values[i] = (i % 3) == 0;
    }

    UIntMask lhs = UIntMask::load_unaligned(uint_values.data());
    FloatMask rhs = FloatMask::load_unaligned(float_values.data());

    auto result = lhs | rhs;

    static_assert(std::same_as<decltype(result), UIntMask>);

    const auto actual = result.to_array();

    for (size_t i = 0; i < UIntMask::size(); ++i) {
        EXPECT_EQ(actual[i], uint_values[i] || float_values[i]);
    }
}


TEST(SimdWidthInteroperabilityTest, DoubleXorUint64MasksCanBeCombined) {
    using DoublePacked = april::simd::Packed<double, 2>;
    using UIntPacked = april::simd::Packed<uint64_t, 2>;

    using DoubleMask = PackedMaskOf<DoublePacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    std::array<bool, DoubleMask::size()> lhs_values{};
    std::array<bool, UIntMask::size()> rhs_values{};

    for (size_t i = 0; i < DoubleMask::size(); ++i) {
        lhs_values[i] = (i % 2) == 0;
        rhs_values[i] = (i % 3) == 0;
    }

    DoubleMask lhs = DoubleMask::load_unaligned(lhs_values.data());
    UIntMask rhs = UIntMask::load_unaligned(rhs_values.data());

    auto result = lhs ^ rhs;

    static_assert(std::same_as<decltype(result), DoubleMask>);

    const auto actual = result.to_array();

    for (size_t i = 0; i < DoubleMask::size(); ++i) {
        EXPECT_EQ(actual[i], lhs_values[i] != rhs_values[i]);
    }
}


// ---------------------------------------------------------
// CONVERTED MASK WITH SELECT
// ---------------------------------------------------------

TEST(SimdWidthInteroperabilityTest, Uint32MaskCanSelectFloat4PackedValues) {
    using FloatPacked = april::simd::Packed<float, 4>;
    using UIntPacked = april::simd::Packed<uint32_t, 4>;

    using FloatMask = PackedMaskOf<FloatPacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    std::array<bool, UIntMask::size()> mask_values{};

    for (size_t i = 0; i < UIntMask::size(); ++i) {
        mask_values[i] = (i % 2) == 0;
    }

    UIntMask source_mask = UIntMask::load_unaligned(mask_values.data());
    FloatMask float_mask = static_cast<FloatMask>(source_mask);

    FloatPacked true_value(10.0f);
    FloatPacked false_value(-10.0f);

    FloatPacked result = select(float_mask, true_value, false_value);
    const auto actual = result.to_array();

    for (size_t i = 0; i < FloatPacked::size(); ++i) {
        EXPECT_EQ(actual[i], mask_values[i] ? 10.0f : -10.0f);
    }
}


TEST(SimdWidthInteroperabilityTest, DoubleMaskCanSelectUint64PackedValues) {
    using DoublePacked = april::simd::Packed<double, 2>;
    using UIntPacked = april::simd::Packed<uint64_t, 2>;

    using DoubleMask = PackedMaskOf<DoublePacked>;
    using UIntMask = PackedMaskOf<UIntPacked>;

    std::array<bool, DoubleMask::size()> mask_values{};

    for (size_t i = 0; i < DoubleMask::size(); ++i) {
        mask_values[i] = (i % 2) == 1;
    }

    DoubleMask source_mask = DoubleMask::load_unaligned(mask_values.data());
    UIntMask uint_mask = static_cast<UIntMask>(source_mask);

    UIntPacked true_value(uint64_t{42});
    UIntPacked false_value(uint64_t{7});

    UIntPacked result = select(uint_mask, true_value, false_value);
    const auto actual = result.to_array();

    for (size_t i = 0; i < UIntPacked::size(); ++i) {
        EXPECT_EQ(actual[i], mask_values[i] ? uint64_t{42} : uint64_t{7});
    }
}


// ---------------------------------------------------------
// DIFFERENT WIDTHS MUST NOT INTEROPERATE
// ---------------------------------------------------------

TEST(SimdWidthInteroperabilityTest, DifferentWidthMasksAreNotConvertible) {
    using Float4 = april::simd::Packed<float, 4>;
    using Double2 = april::simd::Packed<double, 2>;
    using UInt4 = april::simd::Packed<uint32_t, 4>;
    using UInt64x2 = april::simd::Packed<uint64_t, 2>;

    using Float4Mask = PackedMaskOf<Float4>;
    using Double2Mask = PackedMaskOf<Double2>;
    using UInt4Mask = PackedMaskOf<UInt4>;
    using UInt64x2Mask = PackedMaskOf<UInt64x2>;

    static_assert(Float4Mask::size() == 4);
    static_assert(Double2Mask::size() == 2);
    static_assert(UInt4Mask::size() == 4);
    static_assert(UInt64x2Mask::size() == 2);

    static_assert(!std::is_convertible_v<Float4Mask, Double2Mask>);
    static_assert(!std::is_convertible_v<Double2Mask, Float4Mask>);
    static_assert(!std::is_convertible_v<UInt4Mask, UInt64x2Mask>);
    static_assert(!std::is_convertible_v<UInt64x2Mask, UInt4Mask>);

    SUCCEED();
}


TEST(SimdWidthInteroperabilityTest, DifferentWidthPackedTypesRemainDistinct) {
    using Float4 = april::simd::Packed<float, 4>;
    using FloatNative = april::simd::Packed<float>;
    using Double2 = april::simd::Packed<double, 2>;

    static_assert(!std::same_as<Float4, FloatNative>);
    static_assert(!std::same_as<Float4, Double2>);

    static_assert(Float4::size() == 4);
    static_assert(Double2::size() == 2);

    SUCCEED();
}



TYPED_TEST(SimdExplicitWidthTest, MaskRotationUsesExplicitWidth) {
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Width> input{};
    input[0] = true;

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.rotate_right();

    const auto output = mask.to_array();

    for (size_t i = 0; i < TestFixture::Width; ++i) {
        const bool expected = i == 1 % TestFixture::Width;
        EXPECT_EQ(output[i], expected) << "Mismatch at lane " << i;
    }
}





