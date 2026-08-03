#pragma once


namespace april::simd::internal {
    template<typename... Ts>
   struct packed_argument;

    template<typename T, typename... Ts>
    struct packed_argument<T, Ts...> : std::conditional_t<
        IsSimdType<std::remove_cvref_t<T>>,
        std::type_identity<std::remove_cvref_t<T>>,
        packed_argument<Ts...>
    > {};

    template<>
    struct packed_argument<> {};

    template<typename... Ts>
    using packed_argument_t = packed_argument<Ts...>::type;

    template<typename... Ts>
    concept CompatiblePackedArguments = requires {
        typename packed_argument<Ts...>::type;

        requires IsSimdType<packed_argument_t<Ts...>>;
        requires (std::convertible_to<const Ts&, packed_argument_t<Ts...>> && ...);
    };
}

namespace april::simd {
    // Selection
    template<typename M, typename A, typename B>
    requires
        internal::CompatiblePackedArguments<A, B> &&
        internal::IsSimdMaskConvertibleTo<M, typename internal::packed_argument_t<A, B>::mask_type>
    [[nodiscard]] auto select(const M& mask, const A& true_value, const B& false_value) {
        using Packed = internal::packed_argument_t<A, B>;
        using Mask = Packed::mask_type;

        return Packed::select(
            static_cast<Mask>(mask),
            static_cast<Packed>(true_value),
            static_cast<Packed>(false_value)
        );
    }


    // Basic
    template<typename T, size_t W> [[nodiscard]] auto abs(const Packed<T, W>& x) { return Packed<T, W>::abs(x); }

    template<typename A, typename B>
    requires internal::CompatiblePackedArguments<A, B>
    [[nodiscard]] auto min(const A& a, const B& b) {
        using Packed = internal::packed_argument_t<A, B>;
        return Packed::min(static_cast<Packed>(a), static_cast<Packed>(b));
    }

    template<typename A, typename B>
    requires internal::CompatiblePackedArguments<A, B>
    [[nodiscard]] auto max(const A& a, const B& b) {
        using Packed = internal::packed_argument_t<A, B>;
        return Packed::max(static_cast<Packed>(a), static_cast<Packed>(b));
    }

    template<typename X, typename L, typename H>
    requires internal::CompatiblePackedArguments<X, L, H>
    [[nodiscard]] auto clamp(const X& x, const L& lo, const H& hi) {
        using Packed = internal::packed_argument_t<X, L, H>;

        return Packed::clamp(static_cast<Packed>(x), static_cast<Packed>(lo), static_cast<Packed>(hi));
    }


    // Roots and powers
    template<typename T, size_t W> [[nodiscard]] auto sqrt(const Packed<T, W>& x) { return Packed<T, W>::sqrt(x); }
    template<typename T, size_t W> [[nodiscard]] auto rsqrt(const Packed<T, W>& x) { return Packed<T, W>::rsqrt(x); }
    template<typename T, size_t W> [[nodiscard]] auto cbrt(const Packed<T, W>& x) { return Packed<T, W>::cbrt(x); }
    template<typename T, size_t W> [[nodiscard]] auto hypot(const Packed<T, W>& x, const Packed<T, W>& y) { return Packed<T, W>::hypot(x, y); }
    template<typename T, size_t W> [[nodiscard]] auto pow(const Packed<T, W>& x, const Packed<T, W>& y) { return Packed<T, W>::pow(x, y); }

    // Exponential and logarithmic
    template<typename T, size_t W> [[nodiscard]] auto exp(const Packed<T, W>& x) { return Packed<T, W>::exp(x); }
    template<typename T, size_t W> [[nodiscard]] auto exp2(const Packed<T, W>& x) { return Packed<T, W>::exp2(x); }
    template<typename T, size_t W> [[nodiscard]] auto expm1(const Packed<T, W>& x) { return Packed<T, W>::expm1(x); }
    template<typename T, size_t W> [[nodiscard]] auto log(const Packed<T, W>& x) { return Packed<T, W>::log(x); }
    template<typename T, size_t W> [[nodiscard]] auto ln(const Packed<T, W>& x) { return Packed<T, W>::log(x); }
    template<typename T, size_t W> [[nodiscard]] auto log2(const Packed<T, W>& x) { return Packed<T, W>::log2(x); }
    template<typename T, size_t W> [[nodiscard]] auto log10(const Packed<T, W>& x) { return Packed<T, W>::log10(x); }
    template<typename T, size_t W> [[nodiscard]] auto log1p(const Packed<T, W>& x) { return Packed<T, W>::log1p(x); }

    // Trigonometric
    template<typename T, size_t W> [[nodiscard]] auto sin(const Packed<T, W>& x) { return Packed<T, W>::sin(x); }
    template<typename T, size_t W> [[nodiscard]] auto cos(const Packed<T, W>& x) { return Packed<T, W>::cos(x); }
    template<typename T, size_t W> [[nodiscard]] auto sincos(const Packed<T, W>& x) { return Packed<T, W>::sincos(x); }
    template<typename T, size_t W> [[nodiscard]] auto tan(const Packed<T, W>& x) { return Packed<T, W>::tan(x); }
    template<typename T, size_t W> [[nodiscard]] auto asin(const Packed<T, W>& x) { return Packed<T, W>::asin(x); }
    template<typename T, size_t W> [[nodiscard]] auto acos(const Packed<T, W>& x) { return Packed<T, W>::acos(x); }
    template<typename T, size_t W> [[nodiscard]] auto atan(const Packed<T, W>& x) { return Packed<T, W>::atan(x); }

    template<typename Y, typename X>
    requires internal::CompatiblePackedArguments<Y, X>
    [[nodiscard]] auto atan2(const Y& y, const X& x) {
        using Packed = internal::packed_argument_t<Y, X>;
        return Packed::atan2(static_cast<Packed>(y), static_cast<Packed>(x));
    }
    // Hyperbolic
    template<typename T, size_t W> [[nodiscard]] auto sinh(const Packed<T, W>& x) { return Packed<T, W>::sinh(x); }
    template<typename T, size_t W> [[nodiscard]] auto cosh(const Packed<T, W>& x) { return Packed<T, W>::cosh(x); }
    template<typename T, size_t W> [[nodiscard]] auto tanh(const Packed<T, W>& x) { return Packed<T, W>::tanh(x); }
    template<typename T, size_t W> [[nodiscard]] auto asinh(const Packed<T, W>& x) { return Packed<T, W>::asinh(x); }
    template<typename T, size_t W> [[nodiscard]] auto acosh(const Packed<T, W>& x) { return Packed<T, W>::acosh(x); }
    template<typename T, size_t W> [[nodiscard]] auto atanh(const Packed<T, W>& x) { return Packed<T, W>::atanh(x); }

    // Rounding
    template<typename T, size_t W> [[nodiscard]] auto floor(const Packed<T, W>& x) { return Packed<T, W>::floor(x); }
    template<typename T, size_t W> [[nodiscard]] auto ceil(const Packed<T, W>& x) { return Packed<T, W>::ceil(x); }
    template<typename T, size_t W> [[nodiscard]] auto round(const Packed<T, W>& x) { return Packed<T, W>::round(x); }
    template<typename T, size_t W> [[nodiscard]] auto trunc(const Packed<T, W>& x) { return Packed<T, W>::trunc(x); }
    template<typename T, size_t W> [[nodiscard]] auto nearbyint(const Packed<T, W>& x) { return Packed<T, W>::nearbyint(x); }

    // Numeric
    template<typename X, typename Y, typename Z>
    requires internal::CompatiblePackedArguments<X, Y, Z>
    [[nodiscard]] auto fma(const X& x, const Y& y, const Z& z) {
        using Packed = internal::packed_argument_t<X, Y, Z>;

        return Packed::fma(
            static_cast<Packed>(x),
            static_cast<Packed>(y),
            static_cast<Packed>(z)
        );
    }

    template<typename X, typename Y>
    requires internal::CompatiblePackedArguments<X, Y>
    [[nodiscard]] auto fmod(const X& x, const Y& y) {
        using Packed = internal::packed_argument_t<X, Y>;
        return Packed::fmod(static_cast<Packed>(x), static_cast<Packed>(y));
    }

    template<typename X, typename Y>
    requires internal::CompatiblePackedArguments<X, Y>
    [[nodiscard]] auto remainder(const X& x, const Y& y) {
        using Packed = internal::packed_argument_t<X, Y>;
        return Packed::remainder(static_cast<Packed>(x), static_cast<Packed>(y));
    }

    template<typename X, typename Y>
    requires internal::CompatiblePackedArguments<X, Y>
    [[nodiscard]] auto copysign(const X& x, const Y& y) {
        using Packed = internal::packed_argument_t<X, Y>;
        return Packed::copysign(static_cast<Packed>(x), static_cast<Packed>(y));
    }

    // Classification
    template<typename T, size_t W> [[nodiscard]] auto isnan(const Packed<T, W>& x) { return Packed<T, W>::isnan(x); }
    template<typename T, size_t W> [[nodiscard]] auto isinf(const Packed<T, W>& x) { return Packed<T, W>::isinf(x); }
    template<typename T, size_t W> [[nodiscard]] auto isfinite(const Packed<T, W>& x) { return Packed<T, W>::isfinite(x); }
    template<typename T, size_t W> [[nodiscard]] auto signbit(const Packed<T, W>& x) { return Packed<T, W>::signbit(x); }
}

namespace april {
    using simd::select;

    using simd::abs;
    using simd::min;
    using simd::max;
    using simd::clamp;

    using simd::sqrt;
    using simd::rsqrt;
    using simd::cbrt;
    using simd::hypot;
    using simd::pow;

    using simd::exp;
    using simd::exp2;
    using simd::expm1;
    using simd::log;
    using simd::ln;
    using simd::log2;
    using simd::log10;
    using simd::log1p;

    using simd::sin;
    using simd::cos;
    using simd::sincos;
    using simd::tan;
    using simd::asin;
    using simd::acos;
    using simd::atan;
    using simd::atan2;

    using simd::sinh;
    using simd::cosh;
    using simd::tanh;
    using simd::asinh;
    using simd::acosh;
    using simd::atanh;

    using simd::floor;
    using simd::ceil;
    using simd::round;
    using simd::trunc;
    using simd::nearbyint;

    using simd::fma;
    using simd::fmod;
    using simd::remainder;
    using simd::copysign;

    using simd::isnan;
    using simd::isinf;
    using simd::isfinite;
    using simd::signbit;
}