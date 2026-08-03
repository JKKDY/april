#pragma once
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>

#include <april/base/concepts.hpp>

namespace april {
    // Selection
    template<typename T>
    [[nodiscard]] constexpr T select(bool mask, T true_value, T false_value) {
        return mask ? true_value : false_value;
    }

    // Basic
    template<IsScalar T>
    [[nodiscard]] constexpr auto abs(T x) {
        if constexpr (std::is_unsigned_v<T>) return x;
        else return std::abs(x);
    }

    template<IsScalar A, IsScalar B>
    [[nodiscard]] constexpr auto min(A a, B b) {
        using R = std::common_type_t<A, B>;
        return std::min(static_cast<R>(a), static_cast<R>(b));
    }

    template<IsScalar A, IsScalar B>
    [[nodiscard]] constexpr auto max(A a, B b) {
        using R = std::common_type_t<A, B>;
        return std::max(static_cast<R>(a), static_cast<R>(b));
    }

    template<IsScalar T, IsScalar L, IsScalar H>
    [[nodiscard]] constexpr auto clamp(T x, L lo, H hi) {
        using R = std::common_type_t<T, L, H>;
        return std::clamp(static_cast<R>(x), static_cast<R>(lo), static_cast<R>(hi));
    }

    // Roots and powers
    template<IsScalar T> [[nodiscard]] auto sqrt(T x) { return std::sqrt(x); }

    template<IsScalar T>
    [[nodiscard]] auto rsqrt(T x) {
        auto root = std::sqrt(x);
        return decltype(root){1} / root;
    }

    template<IsScalar T> [[nodiscard]] auto cbrt(T x) { return std::cbrt(x); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto hypot(A x, B y) { return std::hypot(x, y); }
    template<IsScalar A, IsScalar B, IsScalar C> [[nodiscard]] auto hypot(A x, B y, C z) { return std::hypot(x, y, z); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto pow(A x, B y) { return std::pow(x, y); }

    // Exponential and logarithmic
    template<IsScalar T> [[nodiscard]] auto exp(T x) { return std::exp(x); }
    template<IsScalar T> [[nodiscard]] auto exp2(T x) { return std::exp2(x); }
    template<IsScalar T> [[nodiscard]] auto expm1(T x) { return std::expm1(x); }
    template<IsScalar T> [[nodiscard]] auto log(T x) { return std::log(x); }
    template<IsScalar T> [[nodiscard]] auto ln(T x) { return log(x); }
    template<IsScalar T> [[nodiscard]] auto log2(T x) { return std::log2(x); }
    template<IsScalar T> [[nodiscard]] auto log10(T x) { return std::log10(x); }
    template<IsScalar T> [[nodiscard]] auto log1p(T x) { return std::log1p(x); }

    // Trigonometric
    template<IsScalar T> [[nodiscard]] auto sin(T x) { return std::sin(x); }
    template<IsScalar T> [[nodiscard]] auto cos(T x) { return std::cos(x); }

    template<IsScalar T>
    [[nodiscard]] auto sincos(T x) {
        auto s = std::sin(x);
        return std::pair{s, std::cos(x)};
    }

    template<IsScalar T> [[nodiscard]] auto tan(T x) { return std::tan(x); }
    template<IsScalar T> [[nodiscard]] auto asin(T x) { return std::asin(x); }
    template<IsScalar T> [[nodiscard]] auto acos(T x) { return std::acos(x); }
    template<IsScalar T> [[nodiscard]] auto atan(T x) { return std::atan(x); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto atan2(A y, B x) { return std::atan2(y, x); }

    // Hyperbolic
    template<IsScalar T> [[nodiscard]] auto sinh(T x) { return std::sinh(x); }
    template<IsScalar T> [[nodiscard]] auto cosh(T x) { return std::cosh(x); }
    template<IsScalar T> [[nodiscard]] auto tanh(T x) { return std::tanh(x); }
    template<IsScalar T> [[nodiscard]] auto asinh(T x) { return std::asinh(x); }
    template<IsScalar T> [[nodiscard]] auto acosh(T x) { return std::acosh(x); }
    template<IsScalar T> [[nodiscard]] auto atanh(T x) { return std::atanh(x); }

    // Rounding
    template<IsScalar T> [[nodiscard]] auto floor(T x) { return std::floor(x); }
    template<IsScalar T> [[nodiscard]] auto ceil(T x) { return std::ceil(x); }
    template<IsScalar T> [[nodiscard]] auto round(T x) { return std::round(x); }
    template<IsScalar T> [[nodiscard]] auto trunc(T x) { return std::trunc(x); }
    template<IsScalar T> [[nodiscard]] auto nearbyint(T x) { return std::nearbyint(x); }

    // Numeric
    template<IsScalar A, IsScalar B, IsScalar C> [[nodiscard]] auto fma(A x, B y, C z) { return std::fma(x, y, z); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto fmod(A x, B y) { return std::fmod(x, y); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto remainder(A x, B y) { return std::remainder(x, y); }
    template<IsScalar A, IsScalar B> [[nodiscard]] auto copysign(A magnitude, B sign) { return std::copysign(magnitude, sign); }

    // Classification
    template<IsScalar T> [[nodiscard]] bool isnan(T x) { return std::isnan(x); }
    template<IsScalar T> [[nodiscard]] bool isinf(T x) { return std::isinf(x); }
    template<IsScalar T> [[nodiscard]] bool isfinite(T x) { return std::isfinite(x); }
    template<IsScalar T> [[nodiscard]] bool signbit(T x) { return std::signbit(x); }
}