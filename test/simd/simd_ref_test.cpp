#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "april/simd/packed_ref.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/backends/backend_scalar.hpp"




using BackendTypes = testing::Types<
    april::simd::Packed<double>,
    april::simd::Packed<float>,
    april::simd::internal::scalar::Packed<double>,
    april::simd::internal::scalar::Packed<float>
>;

template<typename T>
class SimdRefTest : public testing::Test {
public:
    using Packed = T;
    using Scalar = typename T::value_type;
    using Location = april::simd::ContiguousLocation<Scalar>;
    using Ref = april::simd::PackedRef<Location>;

    static_assert(std::same_as<
        typename Location::packed_type,
        Packed
    >);

    std::vector<Scalar> buffer;

    void SetUp() override {
        const size_t safe_size = std::max<size_t>(Packed::size(), 16);
        buffer.resize(safe_size);

        std::fill(
            buffer.begin(),
            buffer.end(),
            static_cast<Scalar>(0.0)
        );
    }

    void ExpectAll(const Packed& w, Scalar expected) {
        const auto arr = w.to_array();
        const double epsilon =
            std::is_same_v<Scalar, float> ? 1e-5 : 1e-12;

        for (const auto v : arr) {
            EXPECT_NEAR(
                static_cast<double>(v),
                static_cast<double>(expected),
                epsilon
            );
        }
    }

    void ExpectMemory(Scalar expected) {
        for (size_t i = 0; i < Packed::size(); ++i) {
            EXPECT_DOUBLE_EQ(buffer[i], expected)
                << "Memory mismatch at index " << i;
        }
    }
};

TYPED_TEST_SUITE(SimdRefTest, BackendTypes);

// Load, Store, Broadcast
TYPED_TEST(SimdRefTest, LoadStoreInteraction) {
    using Packed = TestFixture::Packed;
    using Ref = TestFixture::Ref;
    using Location = TestFixture::Location;

    // 1. Setup Memory: [10, 10, 10, 10...]
    std::fill(this->buffer.begin(), this->buffer.end(), static_cast<TestFixture::Scalar>(10.0));

    // Point Ref to the start of the buffer
    Ref ref(Location{this->buffer.data()});

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

// Mixed Arithmetic (Ref, Wide, Scalar)
TYPED_TEST(SimdRefTest, MixedArithmetic) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;
    using Location = TestFixture::Location;

    // We need TWO independent memory blocks for this test
    std::vector<Scalar> buf_a(Packed::size(), 10.0);
    std::vector<Scalar> buf_b(Packed::size(), 2.0);

    Ref a(Location{buf_a.data()});
    Ref b(Location{buf_b.data()});

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

// Compound Assignments
TYPED_TEST(SimdRefTest, CompoundAssignments) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;
    using Location = TestFixture::Location;

    // Setup Memory: [10, 10...]
    std::fill(this->buffer.begin(), this->buffer.end(), static_cast<TestFixture::Scalar>(10.0));
    Ref r(Location{this->buffer.data()});

    // 1. += Scalar (10 + 2 = 12)
    r += 2.0;
    this->ExpectMemory(12.0);

    // 2. *= Wide (12 * 2 = 24)
    r *= Packed(2.0);
    this->ExpectMemory(24.0);

    // 3. -= Ref (Requires second buffer)
    std::vector<Scalar> buf_other(Packed::size(), 4.0);
    Ref other(Location{buf_other.data()});

    r -= other; // 24 - 4 = 20
    this->ExpectMemory(20.0);
}

// Math Functions (ADL & Implicit Conversion)
TYPED_TEST(SimdRefTest, MathFunctions) {
    using Packed = TestFixture::Packed;
    using Ref = TestFixture::Ref;
    using Location = TestFixture::Location;

    // Setup: [25, 25...]
    std::fill(this->buffer.begin(), this->buffer.end(), static_cast<TestFixture::Scalar>(25.0));
    Ref r(Location{this->buffer.data()});

    // sqrt(Ref) -> calls friend sqrt(SimdRef) -> returns Wide
    Packed root = sqrt(r);
    this->ExpectAll(root, 5.0);

    // min(Ref, Scalar) interaction
    // Relies on Ref converting to Wide, and implicit Wide(Scalar) ctor
    Packed m = april::min(r, 5.0); // min(25, 5) -> 5
    this->ExpectAll(m, 5.0);

    Packed m2 = april::max(r, 5.0); // max(25, 5) -> 25
    this->ExpectAll(m2, 25.0);
}

// Comparisons
TYPED_TEST(SimdRefTest, Comparisons) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Ref = TestFixture::Ref;
    using Location = TestFixture::Location;

    // Buffer A: [10, 10...]
    std::vector<Scalar> buf_a(Packed::size(), 10.0);
    // Buffer B: [20, 20...]
    std::vector<Scalar> buf_b(Packed::size(), 20.0);

    Ref a(Location{buf_a.data()});
    Ref b(Location{buf_b.data()});

    // a < b should return a Mask (all true)
    auto mask = (a < b);

    // Check constraints
    static_assert(!std::is_same_v<decltype(mask), bool>, "Comparison must return SIMD Mask");

    // Verify values (assuming mask behaves like xsimd::batch_bool)
    // We can usually cast mask to boolean for specific lanes or use 'all/any'
    // This depends on your Mask API. Assuming 'all(mask)' exists:
    EXPECT_TRUE(april::all(mask));

    // Scalar comparison: a > 50 (False)
    auto mask2 = (a > 50.0);
    EXPECT_FALSE(april::any(mask2));
}














