#include <gtest/gtest.h>

#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <vector>

#include "april/simd/masked_packed.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/simd/packed_concept.hpp"


using MaskedPackedBackendTypes = testing::Types<
    april::simd::Packed<double>,
    april::simd::Packed<float>
>;


template<typename T>
class MaskedPackedTest : public testing::Test {
public:
    using packed_type = T;
    using scalar_type = packed_type::value_type;
    using mask_type = april::simd::PackedMask<scalar_type>;
    using masked_type = april::simd::MaskedPacked<april::simd::Packed<scalar_type>>;
    using ref_type = april::simd::PackedRef<scalar_type, packed_type>;

    static constexpr std::size_t Size = packed_type::size();

    static constexpr double Epsilon =
        std::is_same_v<scalar_type, float> ? 1e-5 : 1e-12;

    static packed_type MakeSequence(scalar_type start = scalar_type{1}, scalar_type step = scalar_type{1}) {
        std::array<scalar_type, Size> values{};

        for (std::size_t i = 0; i < Size; ++i) {
            values[i] = start + static_cast<scalar_type>(i) * step;
        }

        return packed_type::load_unaligned(values.data());
    }

    static packed_type MakePacked(const std::array<scalar_type, Size>& values) {
        return packed_type::load_unaligned(values.data());
    }

    static mask_type MakeMask(const std::array<bool, Size>& values) {
        return mask_type::load_unaligned(values.data());
    }

    static mask_type AlternatingMask(bool first_lane_active = true) {
        std::array<bool, Size> values{};

        for (std::size_t i = 0; i < Size; ++i) {
            values[i] = ((i % 2) == 0) == first_lane_active;
        }

        return MakeMask(values);
    }

    static mask_type FirstHalfMask() {
        std::array<bool, Size> values{};

        for (std::size_t i = 0; i < Size; ++i) {
            values[i] = i < (Size + 1) / 2;
        }

        return MakeMask(values);
    }

    static std::array<bool, Size> AlternatingMaskArray(bool first_lane_active = true) {
        std::array<bool, Size> values{};

        for (std::size_t i = 0; i < Size; ++i) {
            values[i] = ((i % 2) == 0) == first_lane_active;
        }

        return values;
    }

    static void ExpectPacked(const packed_type& actual, const std::array<scalar_type, Size>& expected) {
        const auto values = actual.to_array();

        for (std::size_t i = 0; i < Size; ++i) {
            EXPECT_NEAR(
                static_cast<double>(values[i]),
                static_cast<double>(expected[i]),
                Epsilon
            ) << "Mismatch at lane " << i;
        }
    }

    static void ExpectAll(const packed_type& actual, scalar_type expected) {
        const auto values = actual.to_array();

        for (std::size_t i = 0; i < Size; ++i) {
            EXPECT_NEAR(
                static_cast<double>(values[i]),
                static_cast<double>(expected),
                Epsilon
            ) << "Mismatch at lane " << i;
        }
    }
};


TYPED_TEST_SUITE(MaskedPackedTest, MaskedPackedBackendTypes);


// ---------------------------------------------------------
// TYPE CONTRACT
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, TypeContract) {
    using PackedT = TestFixture::packed_type;
    using MaskT = TestFixture::mask_type;
    using MaskedT = TestFixture::masked_type;

    static_assert(std::is_constructible_v<MaskedT, const PackedT&, const MaskT&>);
    static_assert(std::is_copy_constructible_v<MaskedT>);
    static_assert(std::is_move_constructible_v<MaskedT>);

    // static_assert(!std::is_copy_assignable_v<MaskedT>);
    // static_assert(!std::is_move_assignable_v<MaskedT>);

    static_assert(std::is_convertible_v<const MaskedT&, PackedT>);
    static_assert(!april::simd::IsSimdType<MaskedT>);
}


// ---------------------------------------------------------
// CONSTRUCTION AND READ ACCESS
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, ConstructFromPacked) {
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto source = TestFixture::MakeSequence();

    MaskedT value(source, mask);
    PackedT result = value;

    TestFixture::ExpectPacked(result, source.to_array());
}


TYPED_TEST(MaskedPackedTest, ConstructFromPackedRef) {
    using Scalar =  TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;
    using RefT =    TestFixture::ref_type;

    std::array<Scalar, TestFixture::Size> memory{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        memory[i] = static_cast<Scalar>(i + 10);
    }

    RefT source(memory.data());
    auto mask = TestFixture::AlternatingMask();

    MaskedT value(source, mask);
    PackedT result = value;

    TestFixture::ExpectPacked(result, memory);
}


TYPED_TEST(MaskedPackedTest, ValueReturnsCompletePackedRegister) {
    using MaskedT = TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto source = TestFixture::MakeSequence();

    MaskedT value(source, mask);

    TestFixture::ExpectPacked(value.value(), source.to_array());
}


// ---------------------------------------------------------
// MASKED ASSIGNMENT
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, PackedAssignmentReplacesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);

    auto initial = TestFixture::MakeSequence(Scalar{1});
    auto replacement = TestFixture::MakeSequence(Scalar{100}, Scalar{10});

    MaskedT value(initial, mask);
    value = replacement;

    const auto initial_values = initial.to_array();
    const auto replacement_values = replacement.to_array();

    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i] ? replacement_values[i] : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, ScalarAssignmentReplacesOnlyActiveLanes) {
    using Scalar = typename TestFixture::scalar_type;
    using MaskedT = typename TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value = Scalar{42};

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i] ? Scalar{42} : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


// ---------------------------------------------------------
// PACKED COMPOUND ASSIGNMENT
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, PackedAdditionUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);

    auto initial = TestFixture::MakeSequence(Scalar{1});
    auto rhs = TestFixture::MakeSequence(Scalar{10});

    MaskedT value(initial, mask);
    value += rhs;

    const auto initial_values = initial.to_array();
    const auto rhs_values = rhs.to_array();

    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] + rhs_values[i]
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, PackedSubtractionUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);

    auto initial = TestFixture::MakeSequence(Scalar{20});
    auto rhs = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value -= rhs;

    const auto initial_values = initial.to_array();
    const auto rhs_values = rhs.to_array();

    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] - rhs_values[i]
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, PackedMultiplicationPreservesInactiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);
    value *= rhs;

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] * Scalar{2}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, PackedDivisionPreservesInactiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{2}, Scalar{2});

    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);
    value /= rhs;

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] / Scalar{2}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


// ---------------------------------------------------------
// SCALAR COMPOUND ASSIGNMENT
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, ScalarAdditionUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value += Scalar{5};

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] + Scalar{5}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, ScalarSubtractionUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{10});

    MaskedT value(initial, mask);
    value -= Scalar{3};

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] - Scalar{3}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, ScalarMultiplicationUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value *= Scalar{3};

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] * Scalar{3}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, ScalarDivisionUpdatesOnlyActiveLanes) {
    using Scalar = TestFixture::scalar_type;
    using MaskedT = TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{4}, Scalar{4});

    MaskedT value(initial, mask);
    value /= Scalar{2};

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = mask_values[i]
            ? initial_values[i] / Scalar{2}
            : initial_values[i];
    }

    TestFixture::ExpectPacked(value, expected);
}


// ---------------------------------------------------------
// ARITHMETIC DROPS THE MASK
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, AdditionReturnsOrdinaryPackedValue) {
    using Scalar = TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto initial = TestFixture::MakeSequence(Scalar{1});
    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);

    auto left_result = value + rhs;
    auto right_result = rhs + value;

    static_assert(std::same_as<decltype(left_result), PackedT>);
    static_assert(std::same_as<decltype(right_result), PackedT>);

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = initial_values[i] + Scalar{2};
    }

    TestFixture::ExpectPacked(left_result, expected);
    TestFixture::ExpectPacked(right_result, expected);
}


TYPED_TEST(MaskedPackedTest, SubtractionReturnsOrdinaryPackedValue) {
    using Scalar = TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto initial = TestFixture::MakeSequence(Scalar{5});
    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);

    auto left_result = value - rhs;
    auto right_result = rhs - value;

    static_assert(std::same_as<decltype(left_result), PackedT>);
    static_assert(std::same_as<decltype(right_result), PackedT>);

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected_left{};
    std::array<Scalar, TestFixture::Size> expected_right{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected_left[i] = initial_values[i] - Scalar{2};
        expected_right[i] = Scalar{2} - initial_values[i];
    }

    TestFixture::ExpectPacked(left_result, expected_left);
    TestFixture::ExpectPacked(right_result, expected_right);
}


TYPED_TEST(MaskedPackedTest, MultiplicationReturnsOrdinaryPackedValue) {
    using Scalar = TestFixture::scalar_type;
    using PackedT = TestFixture::packed_type;
    using MaskedT = TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto initial = TestFixture::MakeSequence(Scalar{1});
    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);

    auto left_result = value * rhs;
    auto right_result = rhs * value;

    static_assert(std::same_as<decltype(left_result), PackedT>);
    static_assert(std::same_as<decltype(right_result), PackedT>);

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = initial_values[i] * Scalar{2};
    }

    TestFixture::ExpectPacked(left_result, expected);
    TestFixture::ExpectPacked(right_result, expected);
}


TYPED_TEST(MaskedPackedTest, DivisionReturnsOrdinaryPackedValue) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;

    auto mask = TestFixture::AlternatingMask();
    auto initial = TestFixture::MakeSequence(Scalar{2}, Scalar{2});
    PackedT rhs(Scalar{2});

    MaskedT value(initial, mask);

    auto left_result = value / rhs;
    auto right_result = rhs / value;

    static_assert(std::same_as<decltype(left_result), PackedT>);
    static_assert(std::same_as<decltype(right_result), PackedT>);

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected_left{};
    std::array<Scalar, TestFixture::Size> expected_right{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected_left[i] = initial_values[i] / Scalar{2};
        expected_right[i] = Scalar{2} / initial_values[i];
    }

    TestFixture::ExpectPacked(left_result, expected_left);
    TestFixture::ExpectPacked(right_result, expected_right);
}


// ---------------------------------------------------------
// EXTERNAL MASK OWNERSHIP
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, ObservesChangesToExternalMask) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;
    using MaskT = typename TestFixture::mask_type;

    std::array<bool, TestFixture::Size> first_mask_values{};
    std::array<bool, TestFixture::Size> second_mask_values{};

    first_mask_values[0] = true;

    if constexpr (TestFixture::Size > 1) {
        second_mask_values[1] = true;
    } else {
        second_mask_values[0] = true;
    }

    MaskT mask = TestFixture::MakeMask(first_mask_values);
    MaskedT value(PackedT(Scalar{0}), mask);

    value += TestFixture::MakeSequence(Scalar{1});

    mask = TestFixture::MakeMask(second_mask_values);
    value += TestFixture::MakeSequence(Scalar{10}, Scalar{10});

    std::array<Scalar, TestFixture::Size> expected{};

    expected[0] = Scalar{1};

    if constexpr (TestFixture::Size > 1) {
        expected[1] = Scalar{20};
    } else {
        expected[0] += Scalar{10};
    }

    TestFixture::ExpectPacked(value, expected);
}


// ---------------------------------------------------------
// ROTATION
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, RotationRotatesDataButNotMask) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;

    std::array<bool, TestFixture::Size> mask_values{};
    mask_values[0] = true;

    auto mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value.rotate_right();
    value += PackedT(Scalar{10});

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        const std::size_t source_index =
            (i + TestFixture::Size - 1) % TestFixture::Size;

        expected[i] = initial_values[source_index];
    }

    expected[0] += Scalar{10};

    TestFixture::ExpectPacked(value, expected);
}


TYPED_TEST(MaskedPackedTest, OwnerCanRotateDataAndReplaceSharedMask) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;
    using MaskT = typename TestFixture::mask_type;

    std::array<bool, TestFixture::Size> mask_values{};
    mask_values[0] = true;

    MaskT mask = TestFixture::MakeMask(mask_values);
    auto initial = TestFixture::MakeSequence(Scalar{1});

    MaskedT value(initial, mask);
    value.rotate_right();

    std::array<bool, TestFixture::Size> rotated_mask_values{};
    rotated_mask_values[(0 + 1) % TestFixture::Size] = true;
    mask = TestFixture::MakeMask(rotated_mask_values);

    value += PackedT(Scalar{10});

    const auto initial_values = initial.to_array();
    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        const std::size_t source_index =
            (i + TestFixture::Size - 1) % TestFixture::Size;

        expected[i] = initial_values[source_index];
    }

    expected[(0 + 1) % TestFixture::Size] += Scalar{10};

    TestFixture::ExpectPacked(value, expected);
}


// ---------------------------------------------------------
// MEMORY COMMIT
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, StoreIntoCommitsRegisterToPackedRef) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;
    using RefT = typename TestFixture::ref_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);

    auto initial = TestFixture::MakeSequence(Scalar{1});
    MaskedT value(initial, mask);

    value += PackedT(Scalar{10});

    std::array<Scalar, TestFixture::Size> destination{};
    RefT destination_ref(destination.data());

    destination_ref = value;

    const auto initial_values = initial.to_array();

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        const Scalar expected = mask_values[i]
            ? initial_values[i] + Scalar{10}
            : initial_values[i];

        EXPECT_NEAR(
            static_cast<double>(destination[i]),
            static_cast<double>(expected),
            TestFixture::Epsilon
        ) << "Mismatch at lane " << i;
    }
}


// ---------------------------------------------------------
// LINEAR ITERATION SEMANTICS
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, LinearIterationUpdatesOnlySelectedParticles) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray();
    auto mask = TestFixture::MakeMask(mask_values);

    auto initial_position = TestFixture::MakeSequence(Scalar{1});
    auto initial_velocity = TestFixture::MakeSequence(Scalar{10}, Scalar{10});

    MaskedT position(initial_position, mask);
    MaskedT velocity(initial_velocity, mask);

    PackedT dt(Scalar{0.5});

    position += static_cast<PackedT>(velocity) * dt;
    velocity += Scalar{2};

    const auto position_values = initial_position.to_array();
    const auto velocity_values = initial_velocity.to_array();

    std::array<Scalar, TestFixture::Size> expected_position{};
    std::array<Scalar, TestFixture::Size> expected_velocity{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected_position[i] = mask_values[i]
            ? position_values[i] + velocity_values[i] * Scalar{0.5}
            : position_values[i];

        expected_velocity[i] = mask_values[i]
            ? velocity_values[i] + Scalar{2}
            : velocity_values[i];
    }

    TestFixture::ExpectPacked(position, expected_position);
    TestFixture::ExpectPacked(velocity, expected_velocity);
}


// ---------------------------------------------------------
// PAIR ITERATION SEMANTICS
// ---------------------------------------------------------

TYPED_TEST(MaskedPackedTest, PairOutputsCanUseIndependentMasks) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;

    std::array<bool, TestFixture::Size> p1_mask_values{};
    std::array<bool, TestFixture::Size> p2_mask_values{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        p1_mask_values[i] = (i % 2) == 1;
        p2_mask_values[i] = i < (TestFixture::Size + 1) / 2;
    }

    auto p1_mask = TestFixture::MakeMask(p1_mask_values);
    auto p2_mask = TestFixture::MakeMask(p2_mask_values);

    MaskedT p1_output(PackedT(Scalar{0}), p1_mask);
    MaskedT p2_output(PackedT(Scalar{0}), p2_mask);

    auto contribution = TestFixture::MakeSequence(Scalar{1});

    p1_output += contribution;
    p2_output -= contribution;

    const auto contribution_values = contribution.to_array();

    std::array<Scalar, TestFixture::Size> expected_p1{};
    std::array<Scalar, TestFixture::Size> expected_p2{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected_p1[i] = p1_mask_values[i]
            ? contribution_values[i]
            : Scalar{0};

        expected_p2[i] = p2_mask_values[i]
            ? -contribution_values[i]
            : Scalar{0};
    }

    TestFixture::ExpectPacked(p1_output, expected_p1);
    TestFixture::ExpectPacked(p2_output, expected_p2);
}


TYPED_TEST(MaskedPackedTest, PairOutputsCanShareOneInteractionMask) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;

    auto mask_values = TestFixture::AlternatingMaskArray(false);
    auto active = TestFixture::MakeMask(mask_values);

    MaskedT p1_output(PackedT(Scalar{0}), active);
    MaskedT p2_output(PackedT(Scalar{0}), active);

    auto lhs = TestFixture::MakeSequence(Scalar{1});
    auto rhs = TestFixture::MakeSequence(
        static_cast<Scalar>(TestFixture::Size),
        Scalar{-1}
    );

    PackedT result = lhs * rhs;

    p1_output += result;
    p2_output -= result;

    const auto result_values = result.to_array();

    std::array<Scalar, TestFixture::Size> expected_p1{};
    std::array<Scalar, TestFixture::Size> expected_p2{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected_p1[i] = mask_values[i]
            ? result_values[i]
            : Scalar{0};

        expected_p2[i] = mask_values[i]
            ? -result_values[i]
            : Scalar{0};
    }

    TestFixture::ExpectPacked(p1_output, expected_p1);
    TestFixture::ExpectPacked(p2_output, expected_p2);
}


TYPED_TEST(MaskedPackedTest, RepeatedPairUpdatesRespectChangingExternalMask) {
    using Scalar = typename TestFixture::scalar_type;
    using PackedT = typename TestFixture::packed_type;
    using MaskedT = typename TestFixture::masked_type;
    using MaskT = typename TestFixture::mask_type;

    std::array<bool, TestFixture::Size> mask_values{};
    mask_values[0] = true;

    MaskT mask = TestFixture::MakeMask(mask_values);
    MaskedT output(PackedT(Scalar{0}), mask);

    auto contribution = TestFixture::MakeSequence(Scalar{1});

    for (std::size_t iteration = 0; iteration < TestFixture::Size; ++iteration) {
        output += contribution;

        std::array<bool, TestFixture::Size> next_mask_values{};
        next_mask_values[(iteration + 1) % TestFixture::Size] = true;
        mask = TestFixture::MakeMask(next_mask_values);
    }

    const auto contribution_values = contribution.to_array();

    std::array<Scalar, TestFixture::Size> expected{};

    for (std::size_t i = 0; i < TestFixture::Size; ++i) {
        expected[i] = contribution_values[i];
    }

    TestFixture::ExpectPacked(output, expected);
}






