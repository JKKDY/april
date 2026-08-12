#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "april/simd/packed_ref.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/backends/backend_scalar.hpp"




using BackendTypes = testing::Types<
    april::simd::Packed<double>,
    april::simd::Packed<float>
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





// ---------------------------------------------------------
// STRIDED LOCATION - STATIC STRIDE
// ---------------------------------------------------------

TYPED_TEST(SimdRefTest, StaticStridedLocationLoad) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;

	static constexpr std::ptrdiff_t Stride = 3 * sizeof(Scalar);
	using Location = april::simd::StridedLocation<Scalar, Stride>;

	static_assert(april::simd::IsLocation<Location>);
	static_assert(april::simd::IsWritableLocation<Location>);
	static_assert(std::same_as<typename Location::packed_type, Packed>);

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});

	for (size_t i = 0; i < Packed::size(); ++i)
		memory[i * 3] = static_cast<Scalar>(i + 1);

	Location location(memory.data());
	const Packed result = location.load();
	const auto values = result.to_array();

	for (size_t i = 0; i < Packed::size(); ++i)
		EXPECT_DOUBLE_EQ(values[i], static_cast<Scalar>(i + 1));
}

TYPED_TEST(SimdRefTest, StaticStridedLocationStore) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;

	static constexpr std::ptrdiff_t Stride = 3 * sizeof(Scalar);
	using Location = april::simd::StridedLocation<Scalar, Stride>;

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});

	std::array<Scalar, Packed::size()> values{};
	for (size_t i = 0; i < Packed::size(); ++i)
		values[i] = static_cast<Scalar>(100 + i);

	const Packed packed = Packed::load(values.data());

	Location location(memory.data());
	location.store(packed);

	for (size_t i = 0; i < Packed::size(); ++i) {
		EXPECT_DOUBLE_EQ(memory[i * 3], values[i]);

		EXPECT_DOUBLE_EQ(memory[i * 3 + 1], Scalar{-1});
		EXPECT_DOUBLE_EQ(memory[i * 3 + 2], Scalar{-1});
	}
}

TYPED_TEST(SimdRefTest, StaticStridedPackedRefLoadStore) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;

	static constexpr std::ptrdiff_t Stride = 3 * sizeof(Scalar);
	using Location = april::simd::StridedLocation<Scalar, Stride>;
	using Ref = april::simd::PackedRef<Location>;

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});

	for (size_t i = 0; i < Packed::size(); ++i)
		memory[i * 3] = Scalar{10};

	Ref ref(Location{memory.data()});

	Packed loaded = ref;
	this->ExpectAll(loaded, Scalar{10});

	ref += Scalar{5};

	for (size_t i = 0; i < Packed::size(); ++i) {
		EXPECT_DOUBLE_EQ(memory[i * 3], Scalar{15});
		EXPECT_DOUBLE_EQ(memory[i * 3 + 1], Scalar{-1});
		EXPECT_DOUBLE_EQ(memory[i * 3 + 2], Scalar{-1});
	}
}


// ---------------------------------------------------------
// STRIDED LOCATION - RUNTIME STRIDE
// ---------------------------------------------------------

TYPED_TEST(SimdRefTest, RuntimeStridedLocationLoad) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::StridedLocation<Scalar>;

	const std::ptrdiff_t stride = 3 * sizeof(Scalar);

	static_assert(april::simd::IsLocation<Location>);
	static_assert(april::simd::IsWritableLocation<Location>);
	static_assert(std::same_as<typename Location::packed_type, Packed>);

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});

	for (size_t i = 0; i < Packed::size(); ++i)
		memory[i * 3] = static_cast<Scalar>(i + 1);

	Location location(memory.data(), stride);
	const Packed result = location.load();
	const auto values = result.to_array();

	for (size_t i = 0; i < Packed::size(); ++i)
		EXPECT_DOUBLE_EQ(values[i], static_cast<Scalar>(i + 1));
}

TYPED_TEST(SimdRefTest, RuntimeStridedLocationStore) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::StridedLocation<Scalar>;

	const std::ptrdiff_t stride = 3 * sizeof(Scalar);

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});

	std::array<Scalar, Packed::size()> values{};
	for (size_t i = 0; i < Packed::size(); ++i)
		values[i] = static_cast<Scalar>(100 + i);

	const Packed packed = Packed::load(values.data());

	Location location(memory.data(), stride);
	location.store(packed);

	for (size_t i = 0; i < Packed::size(); ++i) {
		EXPECT_DOUBLE_EQ(memory[i * 3], values[i]);

		EXPECT_DOUBLE_EQ(memory[i * 3 + 1], Scalar{-1});
		EXPECT_DOUBLE_EQ(memory[i * 3 + 2], Scalar{-1});
	}
}

TYPED_TEST(SimdRefTest, RuntimeStridedPackedRefCompoundAssignment) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::StridedLocation<Scalar>;
	using Ref = april::simd::PackedRef<Location>;

	const std::ptrdiff_t stride = 4 * sizeof(Scalar);

	std::vector<Scalar> memory(Packed::size() * 4, Scalar{-1});

	for (size_t i = 0; i < Packed::size(); ++i)
		memory[i * 4] = Scalar{2};

	Ref ref(Location{memory.data(), stride});

	ref *= Scalar{5};

	for (size_t i = 0; i < Packed::size(); ++i) {
		EXPECT_DOUBLE_EQ(memory[i * 4], Scalar{10});
		EXPECT_DOUBLE_EQ(memory[i * 4 + 1], Scalar{-1});
		EXPECT_DOUBLE_EQ(memory[i * 4 + 2], Scalar{-1});
		EXPECT_DOUBLE_EQ(memory[i * 4 + 3], Scalar{-1});
	}
}


// ---------------------------------------------------------
// GATHER LOCATION
// ---------------------------------------------------------

TYPED_TEST(SimdRefTest, GatherLocationLoad) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::GatherLocation<Scalar>;

	static_assert(april::simd::IsLocation<Location>);
	static_assert(april::simd::IsWritableLocation<Location>);
	static_assert(std::same_as<typename Location::packed_type, Packed>);

	std::vector<Scalar> memory(Packed::size() * 2, Scalar{-1});
	std::array<std::ptrdiff_t, Packed::size()> raw_offsets{};
	std::array<size_t, Packed::size()> indices{};

	for (size_t i = 0; i < Packed::size(); ++i) {
		const size_t index =
			2 * ((i + Packed::size() / 2) % Packed::size());

		indices[i] = index;
		raw_offsets[i] = static_cast<std::ptrdiff_t>(index * sizeof(Scalar));
		memory[index] = static_cast<Scalar>(100 + i);
	}

	const april::simd::ByteOffsets<Packed::size()> offsets(raw_offsets);
	Location location(memory.data(), offsets);

	const Packed result = location.load();
	const auto values = result.to_array();

	for (size_t i = 0; i < Packed::size(); ++i)
		EXPECT_DOUBLE_EQ(values[i], static_cast<Scalar>(100 + i));
}

TYPED_TEST(SimdRefTest, GatherLocationStore) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::GatherLocation<Scalar>;

	std::vector<Scalar> memory(Packed::size() * 2, Scalar{-1});
	std::array<std::ptrdiff_t, Packed::size()> raw_offsets{};
	std::array<size_t, Packed::size()> indices{};
	std::array<Scalar, Packed::size()> values{};

	for (size_t i = 0; i < Packed::size(); ++i) {
		const size_t index =
			2 * ((i + Packed::size() / 2) % Packed::size());

		indices[i] = index;
		raw_offsets[i] = static_cast<std::ptrdiff_t>(index * sizeof(Scalar));
		values[i] = static_cast<Scalar>(200 + i);
	}

	const april::simd::ByteOffsets<Packed::size()> offsets(raw_offsets);
	Location location(memory.data(), offsets);

	location.store(Packed::load(values.data()));

	for (size_t i = 0; i < Packed::size(); ++i)
		EXPECT_DOUBLE_EQ(memory[indices[i]], values[i]);

	for (size_t i = 0; i < memory.size(); ++i) {
		const bool selected =
			std::find(indices.begin(), indices.end(), i) != indices.end();

		if (!selected)
			EXPECT_DOUBLE_EQ(memory[i], Scalar{-1});
	}
}

TYPED_TEST(SimdRefTest, GatherPackedRefLoadStore) {
	using Packed = TestFixture::Packed;
	using Scalar = TestFixture::Scalar;
	using Location = april::simd::GatherLocation<Scalar>;
	using Ref = april::simd::PackedRef<Location>;

	std::vector<Scalar> memory(Packed::size() * 3, Scalar{-1});
	std::array<std::ptrdiff_t, Packed::size()> raw_offsets{};
	std::array<size_t, Packed::size()> indices{};

	for (size_t i = 0; i < Packed::size(); ++i) {
		const size_t index =
			3 * ((i + Packed::size() / 2) % Packed::size()) + 1;

		indices[i] = index;
		raw_offsets[i] = static_cast<std::ptrdiff_t>(index * sizeof(Scalar));
		memory[index] = Scalar{10};
	}

	const april::simd::ByteOffsets<Packed::size()> offsets(raw_offsets);
	Ref ref(Location{memory.data(), offsets});

	Packed loaded = ref;
	this->ExpectAll(loaded, Scalar{10});

	ref += Scalar{7};

	for (size_t i = 0; i < Packed::size(); ++i)
		EXPECT_DOUBLE_EQ(memory[indices[i]], Scalar{17});

	for (size_t i = 0; i < memory.size(); ++i) {
		const bool selected =
			std::find(indices.begin(), indices.end(), i) != indices.end();

		if (!selected)
			EXPECT_DOUBLE_EQ(memory[i], Scalar{-1});
	}
}


// ---------------------------------------------------------
// LOCATION CONSTNESS
// ---------------------------------------------------------

TYPED_TEST(SimdRefTest, LocationConstnessContract) {
	using Scalar = TestFixture::Scalar;

	static constexpr std::ptrdiff_t Stride = 2 * sizeof(Scalar);

	using StaticStrided = april::simd::StridedLocation<Scalar, Stride>;
	using ConstStaticStrided = april::simd::StridedLocation<const Scalar, Stride>;

	using RuntimeStrided = april::simd::StridedLocation<Scalar>;
	using ConstRuntimeStrided = april::simd::StridedLocation<const Scalar>;

	using Gather = april::simd::GatherLocation<Scalar>;
	using ConstGather = april::simd::GatherLocation<const Scalar>;

	static_assert(april::simd::IsWritableLocation<StaticStrided>);
	static_assert(!april::simd::IsWritableLocation<ConstStaticStrided>);

	static_assert(april::simd::IsWritableLocation<RuntimeStrided>);
	static_assert(!april::simd::IsWritableLocation<ConstRuntimeStrided>);

	static_assert(april::simd::IsWritableLocation<Gather>);
	static_assert(!april::simd::IsWritableLocation<ConstGather>);
}








