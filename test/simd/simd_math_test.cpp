#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

#include "april/math/math.hpp"
#include "april/simd/masked_packed.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/simd/backends/backend_scalar.hpp"

#ifndef APRIL_FAST_MATH_ENABLED
    #if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
        #define APRIL_FAST_MATH_ENABLED 1
    #else
        #define APRIL_FAST_MATH_ENABLED 0
    #endif
#endif

namespace {

using FloatingSimdTypes = testing::Types<
    april::simd::Packed<float>,
    april::simd::Packed<double>,
    april::simd::internal::scalar::Packed<double>,
    april::simd::internal::scalar::Packed<float>
>;

template<typename T>
class SimdMathTest : public testing::Test {
public:
    using PackedT = T;
    using Scalar = PackedT::value_type;
    using MaskT = PackedT::mask_type;
    static constexpr std::size_t Size = PackedT::size();
    static constexpr double Tolerance = std::same_as<Scalar, float> ? 5e-5 : 2e-12;
    static constexpr double LooseTolerance = std::same_as<Scalar, float> ? 3e-3 : 2e-11;

    template<typename F>
    static std::array<Scalar, Size> Values(F&& fn) {
        std::array<Scalar, Size> values{};
        for (std::size_t i = 0; i < Size; ++i) values[i] = static_cast<Scalar>(fn(i));
        return values;
    }

    static PackedT Load(const std::array<Scalar, Size>& values) {
        return PackedT::load_unaligned(values.data());
    }

    template<typename F>
    static MaskT MakeMask(F&& fn) {
        std::array<bool, Size> values{};
        for (std::size_t i = 0; i < Size; ++i) values[i] = static_cast<bool>(fn(i));
        return MaskT::load_unaligned(values.data());
    }

    template<typename Expected>
    static void ExpectPacked(
        const PackedT& actual,
        Expected&& expected,
        double tolerance = Tolerance
    ) {
        const auto values = actual.to_array();

        for (std::size_t i = 0; i < Size; ++i) {
            const long double a = static_cast<long double>(values[i]);
            const long double e = static_cast<long double>(expected(i));

#if !APRIL_FAST_MATH_ENABLED
            if (std::isnan(e)) {
                EXPECT_TRUE(std::isnan(a)) << "lane " << i;
                continue;
            }

            if (std::isinf(e)) {
                EXPECT_TRUE(std::isinf(a)) << "lane " << i;
                EXPECT_EQ(std::signbit(a), std::signbit(e)) << "lane " << i;
                continue;
            }
#endif

            const double scale = std::max(1.0, std::abs(static_cast<double>(e)));
            EXPECT_NEAR(
                static_cast<double>(a),
                static_cast<double>(e),
                tolerance * scale
            ) << "lane " << i;
        }
    }

    template<typename Expected>
    static void ExpectMask(const MaskT& actual, Expected&& expected) {
        const auto values = actual.to_array();
        for (std::size_t i = 0; i < Size; ++i) EXPECT_EQ(values[i], static_cast<bool>(expected(i))) << "lane " << i;
    }

    template<typename SimdFn, typename ScalarFn>
    static void CheckUnary(const PackedT& x, const std::array<Scalar, Size>& xv,
                           SimdFn&& simd_fn, ScalarFn&& scalar_fn, double tolerance = Tolerance) {
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(x))>, PackedT>);
        ExpectPacked(simd_fn(x), [&](std::size_t i) { return scalar_fn(xv[i]); }, tolerance);
    }

    template<typename S, typename SimdFn, typename ScalarFn>
    static void CheckBinary(const PackedT& a, const std::array<Scalar, Size>& av,
                            const PackedT& b, const std::array<Scalar, Size>& bv,
                            S scalar_arg, SimdFn&& simd_fn, ScalarFn&& scalar_fn,
                            double tolerance = Tolerance) {
        const Scalar s = static_cast<Scalar>(scalar_arg);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, b))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, scalar_arg))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(scalar_arg, b))>, PackedT>);
        ExpectPacked(simd_fn(a, b), [&](std::size_t i) { return scalar_fn(av[i], bv[i]); }, tolerance);
        ExpectPacked(simd_fn(a, scalar_arg), [&](std::size_t i) { return scalar_fn(av[i], s); }, tolerance);
        ExpectPacked(simd_fn(scalar_arg, b), [&](std::size_t i) { return scalar_fn(s, bv[i]); }, tolerance);
    }

    template<typename SA, typename SB, typename SC, typename SimdFn, typename ScalarFn>
    static void CheckTernary(const PackedT& a, const std::array<Scalar, Size>& av,
                             const PackedT& b, const std::array<Scalar, Size>& bv,
                             const PackedT& c, const std::array<Scalar, Size>& cv,
                             SA scalar_a, SB scalar_b, SC scalar_c,
                             SimdFn&& simd_fn, ScalarFn&& scalar_fn,
                             double tolerance = Tolerance) {
        const Scalar sa = static_cast<Scalar>(scalar_a);
        const Scalar sb = static_cast<Scalar>(scalar_b);
        const Scalar sc = static_cast<Scalar>(scalar_c);

        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, b, c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(scalar_a, b, c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, scalar_b, c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, b, scalar_c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(scalar_a, scalar_b, c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(scalar_a, b, scalar_c))>, PackedT>);
        static_assert(std::same_as<std::remove_cvref_t<decltype(simd_fn(a, scalar_b, scalar_c))>, PackedT>);

        ExpectPacked(simd_fn(a, b, c), [&](std::size_t i) { return scalar_fn(av[i], bv[i], cv[i]); }, tolerance);
        ExpectPacked(simd_fn(scalar_a, b, c), [&](std::size_t i) { return scalar_fn(sa, bv[i], cv[i]); }, tolerance);
        ExpectPacked(simd_fn(a, scalar_b, c), [&](std::size_t i) { return scalar_fn(av[i], sb, cv[i]); }, tolerance);
        ExpectPacked(simd_fn(a, b, scalar_c), [&](std::size_t i) { return scalar_fn(av[i], bv[i], sc); }, tolerance);
        ExpectPacked(simd_fn(scalar_a, scalar_b, c), [&](std::size_t i) { return scalar_fn(sa, sb, cv[i]); }, tolerance);
        ExpectPacked(simd_fn(scalar_a, b, scalar_c), [&](std::size_t i) { return scalar_fn(sa, bv[i], sc); }, tolerance);
        ExpectPacked(simd_fn(a, scalar_b, scalar_c), [&](std::size_t i) { return scalar_fn(av[i], sb, sc); }, tolerance);
    }
};

TYPED_TEST_SUITE(SimdMathTest, FloatingSimdTypes);

TEST(ScalarMaskMathTest, Reductions) {
    EXPECT_TRUE(april::all(true));
    EXPECT_FALSE(april::all(false));
    EXPECT_TRUE(april::any(true));
    EXPECT_FALSE(april::any(false));
    EXPECT_FALSE(april::none(true));
    EXPECT_TRUE(april::none(false));
}

TYPED_TEST(SimdMathTest, MaskReductions) {
    using MaskT = typename TestFixture::MaskT;
    const MaskT all_true(true);
    const MaskT all_false(false);
    const auto alternating = TestFixture::MakeMask([](std::size_t i) { return i % 2 == 0; });

    EXPECT_TRUE(april::all(all_true));
    EXPECT_FALSE(april::any(all_false));
    EXPECT_TRUE(april::none(all_false));
    EXPECT_FALSE(april::all(alternating));
    EXPECT_TRUE(april::any(alternating));
    EXPECT_FALSE(april::none(alternating));
}

TYPED_TEST(SimdMathTest, SelectionAndBasicFunctions) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;
    using MaskT = typename TestFixture::MaskT;

    const auto xv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-4.0, -0.5, 2.0, 9.0}; return v[i % 4];
    });
    const auto yv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-2.0, 0.75, 1.5, 10.0}; return v[i % 4];
    });
    const auto lov = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-3.0, -1.0, 0.0, 2.0}; return v[i % 4];
    });
    const auto hiv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-1.0, 1.0, 4.0, 8.0}; return v[i % 4];
    });
    const PackedT x = TestFixture::Load(xv), y = TestFixture::Load(yv);
    const PackedT lo = TestFixture::Load(lov), hi = TestFixture::Load(hiv);
    const MaskT mask = TestFixture::MakeMask([](std::size_t i) { return i % 2 == 0; });

    static_assert(std::same_as<decltype(april::select(mask, x, y)), PackedT>);
    static_assert(std::same_as<decltype(april::select(mask, x, 2.25)), PackedT>);
    static_assert(std::same_as<decltype(april::select(mask, 2.25, y)), PackedT>);
    TestFixture::ExpectPacked(april::select(mask, x, y), [&](std::size_t i) { return i % 2 == 0 ? xv[i] : yv[i]; });
    TestFixture::ExpectPacked(april::select(mask, x, 2.25), [&](std::size_t i) { return i % 2 == 0 ? xv[i] : Scalar{2.25}; });
    TestFixture::ExpectPacked(april::select(mask, 2.25, y), [&](std::size_t i) { return i % 2 == 0 ? Scalar{2.25} : yv[i]; });

    using CompatibleMask = std::conditional_t<
        std::same_as<Scalar, float>,
        april::simd::PackedMask<int>,
        april::simd::PackedMask<std::size_t>
    >;
    const auto mask_values = mask.to_array();
    const CompatibleMask compatible_mask = CompatibleMask::load_unaligned(mask_values.data());
    TestFixture::ExpectPacked(april::select(compatible_mask, x, y), [&](std::size_t i) { return mask_values[i] ? xv[i] : yv[i]; });

    TestFixture::CheckUnary(x, xv,
        [](const auto& a) { return april::abs(a); },
        [](Scalar a) { return std::abs(a); });

    TestFixture::CheckBinary(x, xv, y, yv, 0.25,
        [](const auto& a, const auto& b) { return april::min(a, b); },
        [](Scalar a, Scalar b) { return std::min(a, b); });

    TestFixture::CheckBinary(x, xv, y, yv, 0.25,
        [](const auto& a, const auto& b) { return april::max(a, b); },
        [](Scalar a, Scalar b) { return std::max(a, b); });

    TestFixture::CheckTernary(x, xv, lo, lov, hi, hiv, 0.5, -1.0, 4.0,
        [](const auto& a, const auto& b, const auto& c) { return april::clamp(a, b, c); },
        [](Scalar a, Scalar b, Scalar c) { return std::clamp(a, b, c); });
}

TYPED_TEST(SimdMathTest, RootsAndPowers) {
    static constexpr double RsqrtTolerance = 1e-3;

    using Scalar = typename TestFixture::Scalar;
    const auto xv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {0.25, 1.0, 2.25, 9.0}; return v[i % 4];
    });
    const auto yv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {0.5, 2.0, 3.0, 4.0}; return v[i % 4];
    });
    const auto ev = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {2.0, 3.0, 0.5, -1.0}; return v[i % 4];
    });
    const auto x = TestFixture::Load(xv), y = TestFixture::Load(yv), exponent = TestFixture::Load(ev);

    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::sqrt(a); }, [](Scalar a) { return std::sqrt(a); });
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::rsqrt(a); }, [](Scalar a) { return Scalar{1} / std::sqrt(a); }, RsqrtTolerance);
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::cbrt(a); }, [](Scalar a) { return std::cbrt(a); });

    TestFixture::CheckBinary(x, xv, y, yv, 2.0,
        [](const auto& a, const auto& b) { return april::hypot(a, b); },
        [](Scalar a, Scalar b) { return std::hypot(a, b); });

    TestFixture::CheckBinary(x, xv, exponent, ev, 2.0,
        [](const auto& a, const auto& b) { return april::pow(a, b); },
        [](Scalar a, Scalar b) { return std::pow(a, b); });
}

TYPED_TEST(SimdMathTest, ExponentialAndLogarithmicFunctions) {
    using Scalar = typename TestFixture::Scalar;
    const auto ev = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-1.0, -0.25, 0.5, 1.25}; return v[i % 4];
    });
    const auto lv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {0.25, 0.5, 2.0, 10.0}; return v[i % 4];
    });
    const auto e = TestFixture::Load(ev), l = TestFixture::Load(lv);

    TestFixture::CheckUnary(e, ev, [](const auto& x) { return april::exp(x); }, [](Scalar x) { return std::exp(x); });
    TestFixture::CheckUnary(e, ev, [](const auto& x) { return april::exp2(x); }, [](Scalar x) { return std::exp2(x); });
    TestFixture::CheckUnary(e, ev, [](const auto& x) { return april::expm1(x); }, [](Scalar x) { return std::expm1(x); });
    TestFixture::CheckUnary(l, lv, [](const auto& x) { return april::log(x); }, [](Scalar x) { return std::log(x); });
    TestFixture::CheckUnary(l, lv, [](const auto& x) { return april::ln(x); }, [](Scalar x) { return std::log(x); });
    TestFixture::CheckUnary(l, lv, [](const auto& x) { return april::log2(x); }, [](Scalar x) { return std::log2(x); });
    TestFixture::CheckUnary(l, lv, [](const auto& x) { return april::log10(x); }, [](Scalar x) { return std::log10(x); });
    TestFixture::CheckUnary(e, ev, [](const auto& x) { return april::log1p(x); }, [](Scalar x) { return std::log1p(x); });
}

TYPED_TEST(SimdMathTest, TrigonometricFunctions) {
    using PackedT = typename TestFixture::PackedT;
    using Scalar = typename TestFixture::Scalar;
    const auto av = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-1.0, -0.5, 0.25, 1.0}; return v[i % 4];
    });
    const auto iv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-0.9, -0.25, 0.5, 0.9}; return v[i % 4];
    });
    const auto yv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-2.0, -1.0, 1.0, 2.0}; return v[i % 4];
    });
    const auto xv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {1.0, -1.0, -2.0, 2.0}; return v[i % 4];
    });
    const auto angle = TestFixture::Load(av), inverse = TestFixture::Load(iv);
    const auto y = TestFixture::Load(yv), x = TestFixture::Load(xv);

    TestFixture::CheckUnary(angle, av, [](const auto& a) { return april::sin(a); }, [](Scalar a) { return std::sin(a); });
    TestFixture::CheckUnary(angle, av, [](const auto& a) { return april::cos(a); }, [](Scalar a) { return std::cos(a); });
    TestFixture::CheckUnary(angle, av, [](const auto& a) { return april::tan(a); }, [](Scalar a) { return std::tan(a); });

    static_assert(std::same_as<decltype(april::sincos(angle)), std::pair<PackedT, PackedT>>);
    const auto [sin_value, cos_value] = april::sincos(angle);
    TestFixture::ExpectPacked(sin_value, [&](std::size_t i) { return std::sin(av[i]); });
    TestFixture::ExpectPacked(cos_value, [&](std::size_t i) { return std::cos(av[i]); });

    TestFixture::CheckUnary(inverse, iv, [](const auto& a) { return april::asin(a); }, [](Scalar a) { return std::asin(a); });
    TestFixture::CheckUnary(inverse, iv, [](const auto& a) { return april::acos(a); }, [](Scalar a) { return std::acos(a); });
    TestFixture::CheckUnary(angle, av, [](const auto& a) { return april::atan(a); }, [](Scalar a) { return std::atan(a); });

    TestFixture::CheckBinary(y, yv, x, xv, 0.75,
        [](const auto& a, const auto& b) { return april::atan2(a, b); },
        [](Scalar a, Scalar b) { return std::atan2(a, b); });
}

TYPED_TEST(SimdMathTest, HyperbolicFunctions) {
    using Scalar = typename TestFixture::Scalar;
    const auto hv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-1.0, -0.25, 0.5, 1.0}; return v[i % 4];
    });
    const auto acv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {1.0, 1.5, 2.0, 4.0}; return v[i % 4];
    });
    const auto atv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-0.75, -0.25, 0.25, 0.75}; return v[i % 4];
    });
    const auto h = TestFixture::Load(hv), ac = TestFixture::Load(acv), at = TestFixture::Load(atv);

    TestFixture::CheckUnary(h, hv, [](const auto& x) { return april::sinh(x); }, [](Scalar x) { return std::sinh(x); });
    TestFixture::CheckUnary(h, hv, [](const auto& x) { return april::cosh(x); }, [](Scalar x) { return std::cosh(x); });
    TestFixture::CheckUnary(h, hv, [](const auto& x) { return april::tanh(x); }, [](Scalar x) { return std::tanh(x); });
    TestFixture::CheckUnary(h, hv, [](const auto& x) { return april::asinh(x); }, [](Scalar x) { return std::asinh(x); });
    TestFixture::CheckUnary(ac, acv, [](const auto& x) { return april::acosh(x); }, [](Scalar x) { return std::acosh(x); });
    TestFixture::CheckUnary(at, atv, [](const auto& x) { return april::atanh(x); }, [](Scalar x) { return std::atanh(x); });
}

TYPED_TEST(SimdMathTest, RoundingFunctions) {
    using Scalar = typename TestFixture::Scalar;
    const auto xv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-2.7, -1.2, 1.4, 2.6}; return v[i % 4];
    });
    const auto x = TestFixture::Load(xv);

    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::floor(a); }, [](Scalar a) { return std::floor(a); });
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::ceil(a); }, [](Scalar a) { return std::ceil(a); });
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::round(a); }, [](Scalar a) { return std::round(a); });
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::trunc(a); }, [](Scalar a) { return std::trunc(a); });
    TestFixture::CheckUnary(x, xv, [](const auto& a) { return april::nearbyint(a); }, [](Scalar a) { return std::nearbyint(a); });
}

TYPED_TEST(SimdMathTest, NumericFunctionsAndMixedArguments) {
    using Scalar = typename TestFixture::Scalar;
    const auto av = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-2.0, -0.5, 1.5, 3.0}; return v[i % 4];
    });
    const auto bv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {0.5, 2.0, -1.0, 4.0}; return v[i % 4];
    });
    const auto cv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {1.0, -2.0, 0.25, 3.0}; return v[i % 4];
    });
    const auto rv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {5.5, -5.5, 7.25, -7.25}; return v[i % 4];
    });
    const auto dv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {2.0, 2.5, -3.0, -2.0}; return v[i % 4];
    });
    const auto sv = TestFixture::Values([](std::size_t i) {
        constexpr double v[] = {-1.0, 1.0, -0.0, 0.0}; return v[i % 4];
    });

    const auto a = TestFixture::Load(av), b = TestFixture::Load(bv), c = TestFixture::Load(cv);
    const auto r = TestFixture::Load(rv), d = TestFixture::Load(dv), signs = TestFixture::Load(sv);

    TestFixture::CheckTernary(a, av, b, bv, c, cv, 1.25, -0.75, 2.5,
        [](const auto& x, const auto& y, const auto& z) { return april::fma(x, y, z); },
        [](Scalar x, Scalar y, Scalar z) { return std::fma(x, y, z); });

    TestFixture::CheckBinary(r, rv, d, dv, 2.0,
        [](const auto& x, const auto& y) { return april::fmod(x, y); },
        [](Scalar x, Scalar y) { return std::fmod(x, y); });

    TestFixture::CheckBinary(r, rv, d, dv, 2.0,
        [](const auto& x, const auto& y) { return april::remainder(x, y); },
        [](Scalar x, Scalar y) { return std::remainder(x, y); });

    TestFixture::CheckBinary(a, av, signs, sv, -1.0,
        [](const auto& x, const auto& y) { return april::copysign(x, y); },
        [](Scalar x, Scalar y) { return std::copysign(x, y); });
}

#if APRIL_FAST_MATH_ENABLED

    TYPED_TEST(SimdMathTest, ClassificationFunctions) {
    GTEST_SKIP() << "NaN and infinity classification requires IEEE floating-point semantics";
}

#else

    TYPED_TEST(SimdMathTest, ClassificationFunctions) {
    using Scalar = typename TestFixture::Scalar;

    const auto values = TestFixture::Values([](std::size_t i) -> Scalar {
        switch (i % 6) {
            case 0: return std::numeric_limits<Scalar>::quiet_NaN();
            case 1: return std::numeric_limits<Scalar>::infinity();
            case 2: return -std::numeric_limits<Scalar>::infinity();
            case 3: return Scalar{-0.0};
            case 4: return Scalar{-2.0};
            default: return Scalar{3.0};
        }
    });

    const auto x = TestFixture::Load(values);

    static_assert(std::same_as<decltype(april::isnan(x)), typename TestFixture::MaskT>);
    static_assert(std::same_as<decltype(april::isinf(x)), typename TestFixture::MaskT>);
    static_assert(std::same_as<decltype(april::isfinite(x)), typename TestFixture::MaskT>);
    static_assert(std::same_as<decltype(april::signbit(x)), typename TestFixture::MaskT>);

    TestFixture::ExpectMask(april::isnan(x), [&](std::size_t i) { return std::isnan(values[i]); });
    TestFixture::ExpectMask(april::isinf(x), [&](std::size_t i) { return std::isinf(values[i]); });
    TestFixture::ExpectMask(april::isfinite(x), [&](std::size_t i) { return std::isfinite(values[i]); });
    TestFixture::ExpectMask(april::signbit(x), [&](std::size_t i) { return std::signbit(values[i]); });
}

#endif

TYPED_TEST(SimdMathTest, PackedRefAndMaskedPackedForwarding) {
	using PackedT = TestFixture::PackedT;
	using Scalar = TestFixture::Scalar;
	using MaskT = TestFixture::MaskT;
	using LocationT = april::simd::ContiguousLocation<Scalar>;
	using RefT = april::simd::PackedRef<LocationT>;
	using MaskedT = april::simd::MaskedPacked<april::simd::Packed<Scalar>>;

	static_assert(std::same_as<typename LocationT::packed_type, PackedT>);

	auto values = TestFixture::Values([](std::size_t i) {
		constexpr double v[] = {0.25, 1.0, 2.25, 4.0}; return v[i % 4];
	});

	RefT ref(LocationT{values.data()});
	const PackedT packed = TestFixture::Load(values);
	const PackedT upper(Scalar{5});
	MaskT mask = TestFixture::MakeMask([](std::size_t i) { return i % 2 == 0; });
	MaskedT masked(packed, mask);

	static_assert(std::same_as<decltype(april::sqrt(ref)), PackedT>);
	static_assert(std::same_as<decltype(april::min(ref, 2.0)), PackedT>);
	static_assert(std::same_as<decltype(april::max(2.0, ref)), PackedT>);
	static_assert(std::same_as<decltype(april::clamp(ref, 0.5, upper)), PackedT>);
	static_assert(std::same_as<decltype(april::fma(ref, 2.0, packed)), PackedT>);
	static_assert(std::same_as<decltype(april::sqrt(masked)), PackedT>);
	static_assert(std::same_as<decltype(april::min(masked, 2.0)), PackedT>);

	TestFixture::ExpectPacked(april::sqrt(ref), [&](std::size_t i) { return std::sqrt(values[i]); });
	TestFixture::ExpectPacked(april::min(ref, 2.0), [&](std::size_t i) { return std::min(values[i], Scalar{2}); });
	TestFixture::ExpectPacked(april::max(2.0, ref), [&](std::size_t i) { return std::max(Scalar{2}, values[i]); });
	TestFixture::ExpectPacked(april::clamp(ref, 0.5, upper), [&](std::size_t i) { return std::clamp(values[i], Scalar{0.5}, Scalar{5}); });
	TestFixture::ExpectPacked(april::fma(ref, 2.0, packed), [&](std::size_t i) { return std::fma(values[i], Scalar{2}, values[i]); });
	TestFixture::ExpectPacked(april::sqrt(masked), [&](std::size_t i) { return std::sqrt(values[i]); });
	TestFixture::ExpectPacked(april::min(masked, 2.0), [&](std::size_t i) { return std::min(values[i], Scalar{2}); });
}

using IntegerSimdTypes = testing::Types<
    april::simd::Packed<int>,
    april::simd::Packed<size_t>,
    april::simd::internal::scalar::Packed<int>,
    april::simd::internal::scalar::Packed<size_t>
>;

template<typename T>
class IntegerSimdMathTest : public testing::Test {
public:
    using PackedT = T;
    using Scalar = PackedT::value_type;
    using MaskT = PackedT::mask_type;
    static constexpr std::size_t Size = PackedT::size();

    template<typename F>
    static std::array<Scalar, Size> Values(F&& fn) {
        std::array<Scalar, Size> values{};
        for (std::size_t i = 0; i < Size; ++i) values[i] = static_cast<Scalar>(fn(i));
        return values;
    }

    static PackedT Load(const std::array<Scalar, Size>& values) { return PackedT::load_unaligned(values.data()); }

    template<typename Expected>
    static void ExpectPacked(const PackedT& actual, Expected&& expected) {
        const auto values = actual.to_array();
        for (std::size_t i = 0; i < Size; ++i) EXPECT_EQ(values[i], static_cast<Scalar>(expected(i))) << "lane " << i;
    }

    template<typename Expected>
    static void ExpectMask(const MaskT& actual, Expected&& expected) {
        const auto values = actual.to_array();
        for (std::size_t i = 0; i < Size; ++i) EXPECT_EQ(values[i], static_cast<bool>(expected(i))) << "lane " << i;
    }
};

TYPED_TEST_SUITE(IntegerSimdMathTest, IntegerSimdTypes);

TYPED_TEST(IntegerSimdMathTest, CommonMathAndScalarInteroperability) {
    using PackedT = TestFixture::PackedT;
    using Scalar = TestFixture::Scalar;
    const auto av = TestFixture::Values([](std::size_t i) {
        if constexpr (std::is_signed_v<Scalar>) {
            constexpr int v[] = {-4, -1, 2, 7}; return v[i % 4];
        } else {
            constexpr std::size_t v[] = {0, 1, 2, 7}; return v[i % 4];
        }
    });
    const auto bv = TestFixture::Values([](std::size_t i) {
        constexpr std::size_t v[] = {3, 2, 5, 1}; return v[i % 4];
    });
    const auto cv = TestFixture::Values([](std::size_t i) {
        constexpr std::size_t v[] = {1, 3, 2, 4}; return v[i % 4];
    });
    const PackedT a = TestFixture::Load(av), b = TestFixture::Load(bv), c = TestFixture::Load(cv);

    TestFixture::ExpectPacked(april::abs(a), [&](std::size_t i) {
        if constexpr (std::is_unsigned_v<Scalar>) return av[i];
        else return static_cast<Scalar>(std::abs(av[i]));
    });

    static_assert(std::same_as<decltype(april::min(a, 3)), PackedT>);
    static_assert(std::same_as<decltype(april::min(3, a)), PackedT>);
    static_assert(std::same_as<decltype(april::max(a, 3)), PackedT>);
    static_assert(std::same_as<decltype(april::clamp(a, 1, b)), PackedT>);
    static_assert(std::same_as<decltype(april::fma(a, 2, c)), PackedT>);

    TestFixture::ExpectPacked(april::min(a, 3), [&](std::size_t i) { return std::min(av[i], Scalar{3}); });
    TestFixture::ExpectPacked(april::min(3, a), [&](std::size_t i) { return std::min(Scalar{3}, av[i]); });
    TestFixture::ExpectPacked(april::max(a, 3), [&](std::size_t i) { return std::max(av[i], Scalar{3}); });
    TestFixture::ExpectPacked(april::max(3, a), [&](std::size_t i) { return std::max(Scalar{3}, av[i]); });
    TestFixture::ExpectPacked(april::clamp(a, 1, b), [&](std::size_t i) { return std::clamp(av[i], Scalar{1}, bv[i]); });
    TestFixture::ExpectPacked(april::fma(a, 2, c), [&](std::size_t i) { return av[i] * Scalar{2} + cv[i]; });
    TestFixture::ExpectPacked(april::fma(2, b, c), [&](std::size_t i) { return Scalar{2} * bv[i] + cv[i]; });

    TestFixture::ExpectPacked(april::floor(a), [&](std::size_t i) { return av[i]; });
    TestFixture::ExpectPacked(april::ceil(a), [&](std::size_t i) { return av[i]; });
    TestFixture::ExpectPacked(april::round(a), [&](std::size_t i) { return av[i]; });
    TestFixture::ExpectPacked(april::trunc(a), [&](std::size_t i) { return av[i]; });
    TestFixture::ExpectPacked(april::nearbyint(a), [&](std::size_t i) { return av[i]; });

    TestFixture::ExpectMask(april::isnan(a), [](std::size_t) { return false; });
    TestFixture::ExpectMask(april::isinf(a), [](std::size_t) { return false; });
    TestFixture::ExpectMask(april::isfinite(a), [](std::size_t) { return true; });
    TestFixture::ExpectMask(april::signbit(a), [&](std::size_t i) {
        if constexpr (std::is_unsigned_v<Scalar>) return false;
        else return av[i] < Scalar{0};
    });
}


TYPED_TEST(SimdMathTest, MaskedPackedInteroperability) {
    using PackedT = TestFixture::PackedT;
    using Scalar = TestFixture::Scalar;
    using MaskT = TestFixture::MaskT;
    using MaskedT = april::simd::MaskedPacked<april::simd::Packed<Scalar>>;

    const auto a_values = TestFixture::Values([](std::size_t i) {
        constexpr double values[] = {0.25, 1.0, 2.25, 4.0};
        return values[i % 4];
    });

    const auto b_values = TestFixture::Values([](std::size_t i) {
        constexpr double values[] = {4.0, 3.0, 2.0, 1.0};
        return values[i % 4];
    });

    const PackedT a = TestFixture::Load(a_values);
    const PackedT b = TestFixture::Load(b_values);
    const MaskT mask = TestFixture::MakeMask([](std::size_t i) {
        return i % 2 == 0;
    });

    const MaskedT masked_a(a, mask);
    const MaskedT masked_b(b, mask);

    static_assert(std::same_as<decltype(april::sqrt(masked_a)), PackedT>);
    static_assert(std::same_as<decltype(april::min(masked_a, b)), PackedT>);
    static_assert(std::same_as<decltype(april::max(a, masked_b)), PackedT>);
    static_assert(std::same_as<decltype(april::clamp(masked_a, 0.5, b)), PackedT>);
    static_assert(std::same_as<decltype(april::fma(masked_a, 2.0, masked_b)), PackedT>);
    static_assert(std::same_as<decltype(april::select(mask, masked_a, b)), PackedT>);

    TestFixture::ExpectPacked(april::sqrt(masked_a), [&](std::size_t i) {
        return std::sqrt(a_values[i]);
    });

    TestFixture::ExpectPacked(april::min(masked_a, b), [&](std::size_t i) {
        return std::min(a_values[i], b_values[i]);
    });

    TestFixture::ExpectPacked(april::max(a, masked_b), [&](std::size_t i) {
        return std::max(a_values[i], b_values[i]);
    });

    TestFixture::ExpectPacked(april::clamp(masked_a, 0.5, b), [&](std::size_t i) {
        return std::clamp(a_values[i], Scalar{0.5}, b_values[i]);
    });

    TestFixture::ExpectPacked(april::fma(masked_a, 2.0, masked_b), [&](std::size_t i) {
        return std::fma(a_values[i], Scalar{2}, b_values[i]);
    });

    TestFixture::ExpectPacked(april::select(mask, masked_a, b), [&](std::size_t i) {
        return i % 2 == 0 ? a_values[i] : b_values[i];
    });
}
}

