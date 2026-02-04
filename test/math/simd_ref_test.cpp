#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "april/simd/backend_std_simd.hpp"
#include "april/simd/backend_xsimd.hpp"
#include "april/simd/packed.hpp"

// Define the Wide types to test
using BackendTypes = testing::Types<
    april::simd::internal::xsimd::Packed<double>,
    april::simd::internal::std_simd::Packed<double>
>;

template <typename T>
class SimdRefTest : public ::testing::Test {
public:
    using Packed = T;
    using Scalar = T::value_type;
    using Ref = april::simd::PackedRef<T>; // Ensure this matches your alias

    // Backing memory must be large enough for at least one SIMD vector
    std::vector<Scalar> buffer;

    void SetUp() override {
        // Resize to safe width (e.g., 8 doubles is enough for AVX-512)
        // We ensure we have enough space for the widest supported register.
        size_t safe_size = std::max<size_t>(Packed::size(), 16);
        buffer.resize(safe_size);

        // Initialize with zeros or a pattern
        std::fill(buffer.begin(), buffer.end(), 0.0);
    }

    // Helper to verify all lanes in a packed match a value
    void ExpectAll(const Packed& w, Scalar expected) {
        auto arr = w.to_array();
        for (auto v : arr) {
            EXPECT_DOUBLE_EQ(v, expected);
        }
    }

    // Helper to verify backing memory was updated correctly
    void ExpectMemory(Scalar expected) {
        for (size_t i = 0; i < Packed::size(); ++i) {
            EXPECT_DOUBLE_EQ(buffer[i], expected)
                << "Memory mismatch at index " << i;
        }
    }
};

TYPED_TEST_SUITE(SimdRefTest, BackendTypes);

// --- 1. Load, Store, Broadcast ---
TYPED_TEST(SimdRefTest, LoadStoreInteraction) {
    using Packed = TestFixture::Packed;
    using Ref = TestFixture::Ref;

    // 1. Setup Memory: [10, 10, 10, 10...]
    std::fill(this->buffer.begin(), this->buffer.end(), 10.0);

    // Point Ref to the start of the buffer
    Ref ref(this->buffer.data());

    // 2. Read (Implicit Load)
    Packed w = ref;
    this->ExpectAll(w, 10.0);

    // 3. Write Scalar (Broadcast & Store)
    // Should write 20.0 to ALL lanes in memory
    ref = 20.0;
    this->ExpectMemory(20.0);

    // 4. Write Wide (Store)
    Packed w2(30.0);
    ref = w2;
    this->ExpectMemory(30.0);
}

// --- 2. Mixed Arithmetic (Ref, Wide, Scalar) ---
TYPED_TEST(SimdRefTest, MixedArithmetic) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;

    // We need TWO independent memory blocks for this test
    std::vector<Scalar> buf_a(Packed::size(), 10.0);
    std::vector<Scalar> buf_b(Packed::size(), 2.0);

    Ref a(buf_a.data());
    Ref b(buf_b.data());

    Packed w(5.0);
    Scalar s = 3.0;

    // 1. Ref + Ref (10 + 2)
    Packed res1 = a + b;
    this->ExpectAll(res1, 12.0);

    // 2. Ref + Scalar (10 + 3)
    Packed res2 = a + s;
    this->ExpectAll(res2, 13.0);

    // 3. Scalar + Ref (3 + 10)
    Packed res3 = s + a;
    this->ExpectAll(res3, 13.0);

    // 4. Ref + Wide (10 + 5)
    Packed res4 = a + w;
    this->ExpectAll(res4, 15.0);

    // 5. Unary Minus (-10)
    Packed res5 = -a;
    this->ExpectAll(res5, -10.0);
}

// --- 3. Compound Assignments ---
TYPED_TEST(SimdRefTest, CompoundAssignments) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;

    // Setup Memory: [10, 10...]
    std::fill(this->buffer.begin(), this->buffer.end(), 10.0);
    Ref r(this->buffer.data());

    // 1. += Scalar (10 + 2 = 12)
    r += 2.0;
    this->ExpectMemory(12.0);

    // 2. *= Wide (12 * 2 = 24)
    r *= Packed(2.0);
    this->ExpectMemory(24.0);

    // 3. -= Ref (Requires second buffer)
    std::vector<Scalar> buf_other(Packed::size(), 4.0);
    Ref other(buf_other.data());

    r -= other; // 24 - 4 = 20
    this->ExpectMemory(20.0);
}

// --- 4. Math Functions (ADL & Implicit Conversion) ---
TYPED_TEST(SimdRefTest, MathFunctions) {
    using Packed = TestFixture::Packed;
    using Ref = TestFixture::Ref;

    // Setup: [25, 25...]
    std::fill(this->buffer.begin(), this->buffer.end(), 25.0);
    Ref r(this->buffer.data());

    // sqrt(Ref) -> calls friend sqrt(SimdRef) -> returns Wide
    Packed root = sqrt(r);
    this->ExpectAll(root, 5.0);

    // min(Ref, Scalar) interaction
    // Relies on Ref converting to Wide, and implicit Wide(Scalar) ctor
    Packed m = min(r, 5.0); // min(25, 5) -> 5
    this->ExpectAll(m, 5.0);

    Packed m2 = max(r, 5.0); // max(25, 5) -> 25
    this->ExpectAll(m2, 25.0);
}

// --- 5. Comparisons (The Mask Check) ---
TYPED_TEST(SimdRefTest, Comparisons) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;

    // Buffer A: [10, 10...]
    std::vector<Scalar> buf_a(Packed::size(), 10.0);
    // Buffer B: [20, 20...]
    std::vector<Scalar> buf_b(Packed::size(), 20.0);

    Ref a(buf_a.data());
    Ref b(buf_b.data());

    // a < b should return a Mask (all true)
    auto mask = (a < b);

    // Check constraints
    static_assert(!std::is_same_v<decltype(mask), bool>, "Comparison must return SIMD Mask");

    // Verify values (assuming mask behaves like xsimd::batch_bool)
    // We can usually cast mask to boolean for specific lanes or use 'all/any'
    // This depends on your Mask API. Assuming 'all(mask)' exists:
    EXPECT_TRUE(all(mask));

    // Scalar comparison: a > 50 (False)
    auto mask2 = (a > 50.0);
    EXPECT_FALSE(any(mask2));
}


