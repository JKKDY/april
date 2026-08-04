#pragma once

#include <concepts>
#include <type_traits>

#include "april/simd/packed_concept.hpp"


namespace april::simd::internal {

    template<typename T>
    concept PackedLike = requires {
        typename std::remove_cvref_t<T>::packed_type;
        requires IsSimdType<typename std::remove_cvref_t<T>::packed_type>;
    };

    // Base case: no packed argument found.
    template<typename... Args>
    struct get_packed {};

    // Skip scalars and unrelated arguments.
    template<typename Head, typename... Tail>
    requires (!PackedLike<Head>)
    struct get_packed<Head, Tail...> : get_packed<Tail...> {};

    // First SIMD value or proxy determines the packed type.
    template<PackedLike Head, typename... Tail>
    struct get_packed<Head, Tail...> {
        using type = std::remove_cvref_t<Head>::packed_type;
    };

    template<typename... Args>
    using get_packed_t = get_packed<Args...>::type;

    template<typename... Args>
    concept CompatiblePackedArguments = requires {
        typename get_packed<Args...>::type;
        requires (
            std::constructible_from<
                typename get_packed<Args...>::type,
                const Args&
            > && ...
        );
    };
}

namespace april::simd {

#define AP_SIMD_UNARY(NAME) \
template<internal::PackedLike X> \
[[nodiscard]] auto NAME(const X& x) { \
using P = internal::get_packed_t<X>; \
return P::NAME(static_cast<P>(x)); \
}

#define AP_SIMD_BINARY(NAME) \
template<typename A, typename B> \
requires internal::CompatiblePackedArguments<A, B> \
[[nodiscard]] auto NAME(const A& a, const B& b) { \
using P = internal::get_packed_t<A, B>; \
return P::NAME(static_cast<P>(a), static_cast<P>(b)); \
}

#define AP_SIMD_TERNARY(NAME) \
template<typename A, typename B, typename C> \
requires internal::CompatiblePackedArguments<A, B, C> \
[[nodiscard]] auto NAME(const A& a, const B& b, const C& c) { \
using P = internal::get_packed_t<A, B, C>; \
return P::NAME(static_cast<P>(a), static_cast<P>(b), static_cast<P>(c)); \
}

    // Selection
    template<typename M, typename A, typename B>
    requires internal::CompatiblePackedArguments<A, B> &&
        internal::IsSimdMaskConvertibleTo<M,typename internal::get_packed_t<A, B>::mask_type>
    [[nodiscard]] auto select(const M& mask, const A& true_value, const B& false_value) {
        using P = internal::get_packed_t<A, B>;
        using Mask = P::mask_type;

        return P::select(
            static_cast<Mask>(mask),
            static_cast<P>(true_value),
            static_cast<P>(false_value)
        );
    }
    AP_SIMD_UNARY(abs)

    AP_SIMD_BINARY(min)
    AP_SIMD_BINARY(max)
    AP_SIMD_TERNARY(clamp)

    AP_SIMD_UNARY(sqrt)
    AP_SIMD_UNARY(rsqrt)
    AP_SIMD_UNARY(cbrt)
    AP_SIMD_BINARY(hypot)
    AP_SIMD_BINARY(pow)

    AP_SIMD_UNARY(exp)
    AP_SIMD_UNARY(exp2)
    AP_SIMD_UNARY(expm1)
    AP_SIMD_UNARY(log)
    AP_SIMD_UNARY(ln)
    AP_SIMD_UNARY(log2)
    AP_SIMD_UNARY(log10)
    AP_SIMD_UNARY(log1p)

    AP_SIMD_UNARY(sin)
    AP_SIMD_UNARY(cos)
    AP_SIMD_UNARY(sincos)
    AP_SIMD_UNARY(tan)
    AP_SIMD_UNARY(asin)
    AP_SIMD_UNARY(acos)
    AP_SIMD_UNARY(atan)
    AP_SIMD_BINARY(atan2)

    AP_SIMD_UNARY(sinh)
    AP_SIMD_UNARY(cosh)
    AP_SIMD_UNARY(tanh)
    AP_SIMD_UNARY(asinh)
    AP_SIMD_UNARY(acosh)
    AP_SIMD_UNARY(atanh)

    AP_SIMD_UNARY(floor)
    AP_SIMD_UNARY(ceil)
    AP_SIMD_UNARY(round)
    AP_SIMD_UNARY(trunc)
    AP_SIMD_UNARY(nearbyint)

    AP_SIMD_TERNARY(fma)
    AP_SIMD_BINARY(fmod)
    AP_SIMD_BINARY(remainder)
    AP_SIMD_BINARY(copysign)

    AP_SIMD_UNARY(isnan)
    AP_SIMD_UNARY(isinf)
    AP_SIMD_UNARY(isfinite)
    AP_SIMD_UNARY(signbit)

    #undef AP_SIMD_TERNARY
    #undef AP_SIMD_BINARY
    #undef AP_SIMD_UNARY

    template<HasMaskReductions M>
    [[nodiscard]] bool all(const M& mask) {
        return std::remove_cvref_t<M>::all(mask);
    }

    template<HasMaskReductions M>
    [[nodiscard]] bool any(const M& mask) {
        return std::remove_cvref_t<M>::any(mask);
    }

    template<HasMaskReductions M>
    [[nodiscard]] bool none(const M& mask) {
        return std::remove_cvref_t<M>::none(mask);
    }
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

    using simd::all;
    using simd::any;
    using simd::none;
}