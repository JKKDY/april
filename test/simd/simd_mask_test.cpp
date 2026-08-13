#include <gtest/gtest.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>

#include "april/simd/packed.hpp"
#include "april/simd/packed_concept.hpp"


namespace {
    using FloatPacked = april::simd::Packed<
        float,
        april::simd::float_width
    >;

    using UInt32Packed = april::simd::Packed<
        uint32_t,
        april::simd::float_width
    >;

    using DoublePacked = april::simd::Packed<
        double,
        april::simd::double_width
    >;

    using UInt64Packed = april::simd::Packed<
        uint64_t,
        april::simd::double_width
    >;

    template<typename PackedT>
    using MaskOf = decltype(
        std::declval<const PackedT&>() ==
        std::declval<const PackedT&>()
    );


    using PublicPackedTypes = testing::Types<
        april::simd::Packed<float, april::simd::float_width>,
        april::simd::Packed<double, april::simd::double_width>,
        april::simd::Packed<uint32_t>,
        april::simd::Packed<uint64_t>
    >;


    template<typename T>
    class SimdMaskTest : public testing::Test {
    public:
        using PackedT = T;
        using Scalar = PackedT::value_type;
        using MaskT = MaskOf<PackedT>;

        static constexpr std::size_t Size = PackedT::size();

        static MaskT MakeAlternatingMask(bool first_active = true) {
            std::array<bool, Size> lanes{};

            for (std::size_t i = 0; i < Size; ++i) {
                lanes[i] = ((i % 2) == 0) == first_active;
            }

            return MaskT::load_unaligned(lanes.data());
        }

        static MaskT MakeSingleLaneMask(std::size_t active_lane) {
            std::array<bool, Size> lanes{};

            if (active_lane < Size) {
                lanes[active_lane] = true;
            }

            return MaskT::load_unaligned(lanes.data());
        }

        static PackedT MakeSequence(Scalar start = Scalar{0}, Scalar step = Scalar{1}) {
            std::array<Scalar, Size> lanes{};

            for (std::size_t i = 0; i < Size; ++i) {
                lanes[i] = start + static_cast<Scalar>(i) * step;
            }

            return PackedT::load_unaligned(lanes.data());
        }

        static void ExpectMask(
            const MaskT& actual,
            const std::array<bool, Size>& expected
        ) {
            const auto lanes = actual.to_array();

            for (std::size_t i = 0; i < Size; ++i) {
                EXPECT_EQ(lanes[i], expected[i])
                    << "Mask mismatch at lane " << i;
            }
        }

        static void ExpectPacked(
            const PackedT& actual,
            const std::array<Scalar, Size>& expected
        ) {
            const auto lanes = actual.to_array();

            for (std::size_t i = 0; i < Size; ++i) {
                EXPECT_EQ(lanes[i], expected[i])
                    << "Packed mismatch at lane " << i;
            }
        }
    };


    template<typename MaskA, typename MaskB>
    void ExpectSameMask(const MaskA& lhs, const MaskB& rhs) {
        static_assert(MaskA::size() == MaskB::size());

        const auto lhs_lanes = lhs.to_array();
        const auto rhs_lanes = rhs.to_array();

        for (std::size_t i = 0; i < MaskA::size(); ++i) {
            EXPECT_EQ(lhs_lanes[i], rhs_lanes[i])
                << "Mask mismatch at lane " << i;
        }
    }


    template<typename SourcePacked, typename TargetPacked>
    void ValidateMaskCast() {
        using SourceMask = MaskOf<SourcePacked>;
        using TargetMask = MaskOf<TargetPacked>;

        if constexpr (SourceMask::size() == TargetMask::size()) {
            static_assert(std::is_convertible_v<SourceMask, TargetMask>);

            std::array<bool, SourceMask::size()> pattern{};

            for (std::size_t i = 0; i < SourceMask::size(); ++i) {
                pattern[i] = (i % 3) != 1;
            }

            SourceMask source = SourceMask::load_unaligned(pattern.data());
            TargetMask target = static_cast<TargetMask>(source);

            ExpectSameMask(source, target);
        } else {
            static_assert(!std::is_convertible_v<SourceMask, TargetMask>);

            SUCCEED()
                << "Mask cast correctly rejected because lane counts differ: "
                << SourceMask::size() << " versus " << TargetMask::size();
        }
    }

}


TYPED_TEST_SUITE(SimdMaskTest, PublicPackedTypes);


// ---------------------------------------------------------
// TYPE CONTRACT
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, SatisfiesPublicSimdContracts) {
    using PackedT = TestFixture::PackedT;
    using MaskT = TestFixture::MaskT;

    static_assert(april::simd::IsSimdType<PackedT>);
    static_assert(april::simd::IsSimdMask<MaskT>);

    static_assert(PackedT::size() == MaskT::size());
    static_assert(std::same_as<decltype(PackedT{} == PackedT{}), MaskT>);
    static_assert(std::same_as<decltype(PackedT{} != PackedT{}), MaskT>);
    static_assert(std::same_as<decltype(PackedT{} < PackedT{}), MaskT>);
    static_assert(std::same_as<decltype(PackedT{} <= PackedT{}), MaskT>);
    static_assert(std::same_as<decltype(PackedT{} > PackedT{}), MaskT>);
    static_assert(std::same_as<decltype(PackedT{} >= PackedT{}), MaskT>);
}


// ---------------------------------------------------------
// BOOLEAN CONSTRUCTION AND REDUCTIONS
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, BooleanConstructorBroadcastsToAllLanes) {
    using MaskT = TestFixture::MaskT;

    MaskT all_true(true);
    MaskT all_false(false);

    EXPECT_TRUE(april::all(all_true));
    EXPECT_TRUE(april::any(all_true));
    EXPECT_FALSE(april::none(all_true));

    EXPECT_FALSE(april::all(all_false));
    EXPECT_FALSE(april::any(all_false));
    EXPECT_TRUE(april::none(all_false));
}


TYPED_TEST(SimdMaskTest, MixedMaskReductionsAreCorrect) {
    using MaskT = TestFixture::MaskT;

    if constexpr (TestFixture::Size > 1) {
        MaskT mask = TestFixture::MakeSingleLaneMask(0);

        EXPECT_FALSE(april::all(mask));
        EXPECT_TRUE(april::any(mask));
        EXPECT_FALSE(april::none(mask));
    } else {
        GTEST_SKIP() << "SIMD mask contains only one lane";
    }
}


// ---------------------------------------------------------
// LOAD, STORE, AND EXPORT
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, LoadStoreRoundTripPreservesLanes) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    std::array<bool, TestFixture::Size> output{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        input[i] = (i % 3) == 0;
    }

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.store_unaligned(output.data());

    EXPECT_EQ(input, output);
}


TYPED_TEST(SimdMaskTest, ToArrayPreservesLanes) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = i < (TestFixture::Size + 1) / 2;
    }

    MaskT mask = MaskT::load_unaligned(expected.data());

    TestFixture::ExpectMask(mask, expected);
}


TYPED_TEST(SimdMaskTest, ToBitmaskMatchesActiveLanes) {
    using MaskT = TestFixture::MaskT;

    static_assert(TestFixture::Size <= 64);

    std::array<bool, TestFixture::Size> lanes{};
    uint64_t expected_bits = 0;

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lanes[i] = (i % 2) == 0;

        if (lanes[i]) {
            expected_bits |= uint64_t{1} << i;
        }
    }

    MaskT mask = MaskT::load_unaligned(lanes.data());

    constexpr uint64_t valid_lane_bits = [] {
        if constexpr (TestFixture::Size == 64) {
            return ~uint64_t{0};
        } else {
            return (uint64_t{1} << TestFixture::Size) - 1;
        }
    }();

    EXPECT_EQ(mask.to_bitmask() & valid_lane_bits, expected_bits);
}


TYPED_TEST(SimdMaskTest, ToStringContainsAllLaneValues) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> lanes{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lanes[i] = (i % 2) == 0;
    }

    MaskT mask = MaskT::load_unaligned(lanes.data());
    const std::string text = mask.to_string();

    EXPECT_FALSE(text.empty());
    EXPECT_EQ(text.front(), '[');
    EXPECT_EQ(text.back(), ']');

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        const std::string expected = lanes[i] ? "true" : "false";
        EXPECT_NE(text.find(expected), std::string::npos);
    }
}


// ---------------------------------------------------------
// SAME-TYPE MASK OPERATORS
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, BitwiseOperatorsAreLaneWise) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> lhs_values{};
    std::array<bool, TestFixture::Size> rhs_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lhs_values[i] = (i % 2) == 0;
        rhs_values[i] = (i % 3) == 0;
    }

    MaskT lhs = MaskT::load_unaligned(lhs_values.data());
    MaskT rhs = MaskT::load_unaligned(rhs_values.data());

    const auto and_values = (lhs & rhs).to_array();
    const auto or_values = (lhs | rhs).to_array();
    const auto xor_values = (lhs ^ rhs).to_array();
    const auto complement_values = (~lhs).to_array();
    const auto logical_not_values = (!lhs).to_array();

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        EXPECT_EQ(and_values[i], lhs_values[i] && rhs_values[i]);
        EXPECT_EQ(or_values[i], lhs_values[i] || rhs_values[i]);
        EXPECT_EQ(xor_values[i], lhs_values[i] != rhs_values[i]);
        EXPECT_EQ(complement_values[i], !lhs_values[i]);
        EXPECT_EQ(logical_not_values[i], !lhs_values[i]);
    }
}


TYPED_TEST(SimdMaskTest, LogicalOperatorsAreLaneWise) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> lhs_values{};
    std::array<bool, TestFixture::Size> rhs_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lhs_values[i] = (i % 2) == 0;
        rhs_values[i] = i < TestFixture::Size / 2;
    }

    MaskT lhs = MaskT::load_unaligned(lhs_values.data());
    MaskT rhs = MaskT::load_unaligned(rhs_values.data());

    const auto and_values = (lhs && rhs).to_array();
    const auto or_values = (lhs || rhs).to_array();

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        EXPECT_EQ(and_values[i], lhs_values[i] && rhs_values[i]);
        EXPECT_EQ(or_values[i], lhs_values[i] || rhs_values[i]);
    }
}


TYPED_TEST(SimdMaskTest, EqualityOperatorsAreLaneWise) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> lhs_values{};
    std::array<bool, TestFixture::Size> rhs_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lhs_values[i] = (i % 2) == 0;
        rhs_values[i] = (i % 3) == 0;
    }

    MaskT lhs = MaskT::load_unaligned(lhs_values.data());
    MaskT rhs = MaskT::load_unaligned(rhs_values.data());

    const auto equal_values = (lhs == rhs).to_array();
    const auto unequal_values = (lhs != rhs).to_array();

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        EXPECT_EQ(equal_values[i], lhs_values[i] == rhs_values[i]);
        EXPECT_EQ(unequal_values[i], lhs_values[i] != rhs_values[i]);
    }
}


// ---------------------------------------------------------
// PACKED COMPARISONS
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, PackedComparisonsProduceExpectedMask) {
    using PackedT = TestFixture::PackedT;
    using Scalar = TestFixture::Scalar;

    std::array<Scalar, TestFixture::Size> lhs_values{};
    std::array<Scalar, TestFixture::Size> rhs_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        lhs_values[i] = static_cast<Scalar>(i);
        rhs_values[i] = static_cast<Scalar>(
            (i % 2) == 0 ? i : i + 1
        );
    }

    PackedT lhs = PackedT::load_unaligned(lhs_values.data());
    PackedT rhs = PackedT::load_unaligned(rhs_values.data());

    const auto equal_values = (lhs == rhs).to_array();
    const auto unequal_values = (lhs != rhs).to_array();
    const auto less_values = (lhs < rhs).to_array();
    const auto less_equal_values = (lhs <= rhs).to_array();
    const auto greater_values = (lhs > rhs).to_array();
    const auto greater_equal_values = (lhs >= rhs).to_array();

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        EXPECT_EQ(equal_values[i], lhs_values[i] == rhs_values[i]);
        EXPECT_EQ(unequal_values[i], lhs_values[i] != rhs_values[i]);
        EXPECT_EQ(less_values[i], lhs_values[i] < rhs_values[i]);
        EXPECT_EQ(less_equal_values[i], lhs_values[i] <= rhs_values[i]);
        EXPECT_EQ(greater_values[i], lhs_values[i] > rhs_values[i]);
        EXPECT_EQ(greater_equal_values[i], lhs_values[i] >= rhs_values[i]);
    }
}


// ---------------------------------------------------------
// SELECT
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, SelectUsesMaskLaneByLane) {
    using PackedT = TestFixture::PackedT;
    using Scalar = TestFixture::Scalar;
    using MaskT = TestFixture::MaskT;

    std::array<Scalar, TestFixture::Size> true_values{};
    std::array<Scalar, TestFixture::Size> false_values{};
    std::array<Scalar, TestFixture::Size> expected{};
    std::array<bool, TestFixture::Size> mask_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        true_values[i] = static_cast<Scalar>(100 + i);
        false_values[i] = static_cast<Scalar>(200 + i);
        mask_values[i] = (i % 2) == 0;

        expected[i] = mask_values[i]
            ? true_values[i]
            : false_values[i];
    }

    PackedT true_value = PackedT::load_unaligned(true_values.data());
    PackedT false_value = PackedT::load_unaligned(false_values.data());
    MaskT mask = MaskT::load_unaligned(mask_values.data());

    PackedT result = april::select(mask, true_value, false_value);

    TestFixture::ExpectPacked(result, expected);
}


// ---------------------------------------------------------
// MASK CASTING CONTRACT
// ---------------------------------------------------------
TEST(SimdMaskCastTest, FloatToUint32PreservesLanes) {
	ValidateMaskCast<FloatPacked, UInt32Packed>();
}

TEST(SimdMaskCastTest, Uint32ToFloatPreservesLanes) {
	ValidateMaskCast<UInt32Packed, FloatPacked>();
}

TEST(SimdMaskCastTest, DoubleToUint64PreservesLanes) {
	ValidateMaskCast<DoublePacked, UInt64Packed>();
}

TEST(SimdMaskCastTest, Uint64ToDoublePreservesLanes) {
	ValidateMaskCast<UInt64Packed, DoublePacked>();
}

TEST(SimdMaskCastTest, FloatToDoubleFollowsLaneCountConstraint) {
	ValidateMaskCast<FloatPacked, DoublePacked>();
}

TEST(SimdMaskCastTest, DoubleToFloatFollowsLaneCountConstraint) {
	ValidateMaskCast<DoublePacked, FloatPacked>();
}

TEST(SimdMaskCastTest, Uint32ToUint64FollowsLaneCountConstraint) {
	ValidateMaskCast<UInt32Packed, UInt64Packed>();
}

TEST(SimdMaskCastTest, Uint64ToUint32FollowsLaneCountConstraint) {
	ValidateMaskCast<UInt64Packed, UInt32Packed>();
}

TEST(SimdMaskCastTest, SameWidthRoundTripPreservesLanes) {
	using FloatMask = MaskOf<FloatPacked>;
	using UIntMask = MaskOf<UInt32Packed>;

	static_assert(FloatMask::size() == UIntMask::size());

	std::array<bool, FloatMask::size()> pattern{};
	for (std::size_t i = 0; i < FloatMask::size(); ++i)
		pattern[i] = (i % 3) == 0;

	FloatMask original = FloatMask::load_unaligned(pattern.data());
	UIntMask converted = static_cast<UIntMask>(original);
	FloatMask round_trip = static_cast<FloatMask>(converted);

	ExpectSameMask(original, converted);
	ExpectSameMask(original, round_trip);
}


// ---------------------------------------------------------
// MIXED-TYPE MASK OPERATORS
// ---------------------------------------------------------

TEST(SimdMixedMaskTest, FloatAndUint32ReturnsLeftMaskType) {
	using FloatMask = MaskOf<FloatPacked>;
	using UIntMask = MaskOf<UInt32Packed>;

	static_assert(FloatMask::size() == UIntMask::size());

	std::array<bool, FloatMask::size()> lhs_values{};
	std::array<bool, UIntMask::size()> rhs_values{};

	for (std::size_t i = 0; i < FloatMask::size(); ++i) {
		lhs_values[i] = (i % 2) == 0;
		rhs_values[i] = (i % 3) == 0;
	}

	FloatMask lhs = FloatMask::load_unaligned(lhs_values.data());
	UIntMask rhs = UIntMask::load_unaligned(rhs_values.data());

	auto result = lhs & rhs;
	static_assert(std::same_as<decltype(result), FloatMask>);

	const auto result_values = result.to_array();
	for (std::size_t i = 0; i < FloatMask::size(); ++i)
		EXPECT_EQ(result_values[i], lhs_values[i] && rhs_values[i]);
}

TEST(SimdMixedMaskTest, Uint32OrFloatReturnsLeftMaskType) {
	using FloatMask = MaskOf<FloatPacked>;
	using UIntMask = MaskOf<UInt32Packed>;

	static_assert(FloatMask::size() == UIntMask::size());

	std::array<bool, FloatMask::size()> lhs_values{};
	std::array<bool, UIntMask::size()> rhs_values{};

	for (std::size_t i = 0; i < UIntMask::size(); ++i) {
		lhs_values[i] = (i % 2) == 0;
		rhs_values[i] = (i % 3) == 0;
	}

	UIntMask lhs = UIntMask::load_unaligned(rhs_values.data());
	FloatMask rhs = FloatMask::load_unaligned(lhs_values.data());

	auto result = lhs | rhs;
	static_assert(std::same_as<decltype(result), UIntMask>);

	const auto result_values = result.to_array();
	for (std::size_t i = 0; i < UIntMask::size(); ++i)
		EXPECT_EQ(result_values[i], rhs_values[i] || lhs_values[i]);
}

TEST(SimdMixedMaskTest, DoubleXorUint64ReturnsLeftMaskType) {
	using DoubleMask = MaskOf<DoublePacked>;
	using UIntMask = MaskOf<UInt64Packed>;

	static_assert(DoubleMask::size() == UIntMask::size());

	std::array<bool, DoubleMask::size()> lhs_values{};
	std::array<bool, UIntMask::size()> rhs_values{};

	for (std::size_t i = 0; i < DoubleMask::size(); ++i) {
		lhs_values[i] = (i % 2) == 0;
		rhs_values[i] = (i % 3) == 0;
	}

	DoubleMask lhs = DoubleMask::load_unaligned(lhs_values.data());
	UIntMask rhs = UIntMask::load_unaligned(rhs_values.data());

	auto result = lhs ^ rhs;
	static_assert(std::same_as<decltype(result), DoubleMask>);

	const auto result_values = result.to_array();
	for (std::size_t i = 0; i < DoubleMask::size(); ++i)
		EXPECT_EQ(result_values[i], lhs_values[i] != rhs_values[i]);
}


// ---------------------------------------------------------
// CONVERTED MASK WITH PACKED SELECT
// ---------------------------------------------------------

TEST(SimdMaskPackedInteractionTest, ConvertedMaskCanSelectFloatLanes) {
	using FloatMask = MaskOf<FloatPacked>;
	using UIntMask = MaskOf<UInt32Packed>;

	static_assert(FloatPacked::size() == UInt32Packed::size());

	std::array<bool, UIntMask::size()> mask_values{};
	for (std::size_t i = 0; i < UIntMask::size(); ++i)
		mask_values[i] = (i % 2) == 0;

	UIntMask source_mask = UIntMask::load_unaligned(mask_values.data());
	FloatMask float_mask = static_cast<FloatMask>(source_mask);

	FloatPacked true_value(1.0f);
	FloatPacked false_value(-1.0f);

	FloatPacked result = april::select(float_mask, true_value, false_value);
	const auto result_values = result.to_array();

	for (std::size_t i = 0; i < FloatPacked::size(); ++i)
		EXPECT_EQ(result_values[i], mask_values[i] ? 1.0f : -1.0f);
}

TEST(SimdMaskPackedInteractionTest, ConvertedMaskCanSelectUint64Lanes) {
	using DoubleMask = MaskOf<DoublePacked>;
	using UIntMask = MaskOf<UInt64Packed>;

	static_assert(DoublePacked::size() == UInt64Packed::size());

	std::array<bool, DoubleMask::size()> mask_values{};
	for (std::size_t i = 0; i < DoubleMask::size(); ++i)
		mask_values[i] = (i % 3) != 0;

	DoubleMask source_mask = DoubleMask::load_unaligned(mask_values.data());
	UIntMask uint_mask = static_cast<UIntMask>(source_mask);

	UInt64Packed true_value(uint64_t{42});
	UInt64Packed false_value(uint64_t{7});

	UInt64Packed result = april::select(uint_mask, true_value, false_value);
	const auto result_values = result.to_array();

	for (std::size_t i = 0; i < UInt64Packed::size(); ++i)
		EXPECT_EQ(result_values[i], mask_values[i] ? uint64_t{42} : uint64_t{7});
}


// ---------------------------------------------------------
// MASK ROTATION
// ---------------------------------------------------------

TYPED_TEST(SimdMaskTest, RotateLeftByOne) {
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    for (size_t i = 0; i < TestFixture::Size; ++i)
        input[i] = (i % 3) == 0;

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.rotate_left();

    const auto output = mask.to_array();

    for (size_t i = 0; i < TestFixture::Size; ++i) {
        const bool expected = input[(i + 1) % TestFixture::Size];
        EXPECT_EQ(output[i], expected) << "Mismatch at lane " << i;
    }
}


TYPED_TEST(SimdMaskTest, RotateRightByOne) {
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    for (size_t i = 0; i < TestFixture::Size; ++i)
        input[i] = (i % 3) == 0;

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.rotate_right();

    const auto output = mask.to_array();

    for (size_t i = 0; i < TestFixture::Size; ++i) {
        const size_t source = (i + TestFixture::Size - 1) % TestFixture::Size;
        EXPECT_EQ(output[i], input[source]) << "Mismatch at lane " << i;
    }
}


TYPED_TEST(SimdMaskTest, RotateLeftByTwo) {
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    for (size_t i = 0; i < TestFixture::Size; ++i)
        input[i] = (i % 3) == 1;

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.template rotate_left<2>();

    const auto output = mask.to_array();

    for (size_t i = 0; i < TestFixture::Size; ++i) {
        const bool expected = input[(i + 2) % TestFixture::Size];
        EXPECT_EQ(output[i], expected) << "Mismatch at lane " << i;
    }
}


TYPED_TEST(SimdMaskTest, RotateRightByTwo) {
    using MaskT = typename TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    for (size_t i = 0; i < TestFixture::Size; ++i)
        input[i] = (i % 3) == 1;

    MaskT mask = MaskT::load_unaligned(input.data());
    mask.template rotate_right<2>();

    const auto output = mask.to_array();

    for (size_t i = 0; i < TestFixture::Size; ++i) {
        const size_t source =
            (i + TestFixture::Size - (2 % TestFixture::Size)) % TestFixture::Size;

        EXPECT_EQ(output[i], input[source]) << "Mismatch at lane " << i;
    }
}


TYPED_TEST(SimdMaskTest, RotateByZeroIsIdentity) {
    using MaskT = typename TestFixture::MaskT;

    MaskT left = TestFixture::MakeAlternatingMask();
    MaskT right = left;

    const auto expected = left.to_array();

    left.template rotate_left<0>();
    right.template rotate_right<0>();

    EXPECT_EQ(left.to_array(), expected);
    EXPECT_EQ(right.to_array(), expected);
}


TYPED_TEST(SimdMaskTest, RotateByFullWidthIsIdentity) {
    using MaskT = typename TestFixture::MaskT;

    MaskT left = TestFixture::MakeAlternatingMask();
    MaskT right = left;

    const auto expected = left.to_array();

    left.template rotate_left<TestFixture::Size>();
    right.template rotate_right<TestFixture::Size>();

    EXPECT_EQ(left.to_array(), expected);
    EXPECT_EQ(right.to_array(), expected);
}


TYPED_TEST(SimdMaskTest, LeftAndRightRotationAreInverses) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> lanes{};
    for (size_t i = 0; i < TestFixture::Size; ++i)
        lanes[i] = (i % 3) != 1;

    MaskT mask = MaskT::load_unaligned(lanes.data());
    const auto expected = mask.to_array();

    mask.template rotate_left<2>();
    mask.template rotate_right<2>();

    EXPECT_EQ(mask.to_array(), expected);
}


TYPED_TEST(SimdMaskTest, RepeatedRightRotationCompletesFullCycle) {
    using MaskT = TestFixture::MaskT;

    std::array<bool, TestFixture::Size> input{};
    input[0] = true;

    MaskT mask = MaskT::load_unaligned(input.data());

    for (size_t i = 0; i < TestFixture::Size; ++i)
        mask.rotate_right();

    EXPECT_EQ(mask.to_array(), input);
}