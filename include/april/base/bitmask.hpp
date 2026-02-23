#pragma once

#include <type_traits>

#define AP_ENABLE_BITMASK_OPERATORS(T)                                         \
    using T##_u = std::underlying_type_t<T>;                                   \
                                                                               \
    /* bitwise logic */                                                        \
    [[nodiscard]] inline constexpr T operator|(T a, T b) noexcept {            \
        return static_cast<T>(static_cast<T##_u>(a) | static_cast<T##_u>(b));  \
    }                                                                          \
    [[nodiscard]] inline constexpr T operator&(T a, T b) noexcept {            \
        return static_cast<T>(static_cast<T##_u>(a) & static_cast<T##_u>(b));  \
    }                                                                          \
    [[nodiscard]] inline constexpr T operator^(T a, T b) noexcept {            \
        return static_cast<T>(static_cast<T##_u>(a) ^ static_cast<T##_u>(b));  \
    }                                                                          \
    [[nodiscard]] inline constexpr T operator~(T a) noexcept {                 \
        return static_cast<T>(~static_cast<T##_u>(a));                         \
    }                                                                          \
                                                                               \
    /* assignment operators */                                                 \
    inline constexpr T& operator|=(T& a, T b) noexcept { return a = a | b; }   \
    inline constexpr T& operator&=(T& a, T b) noexcept { return a = a & b; }   \
    inline constexpr T& operator^=(T& a, T b) noexcept { return a = a ^ b; }   \
                                                                               \
    /* set operations */                                                       \
    [[nodiscard]] inline constexpr T operator-(T a, T b) noexcept {            \
        return a & ~b;                                                         \
    }                                                                          \
    inline constexpr T& operator-=(T& a, T b) noexcept { return a = a & ~b; }  \
                                                                               \
    /* utility checks */                                                       \
    [[nodiscard]] inline constexpr bool operator!(T a) noexcept {              \
        return static_cast<T##_u>(a) == 0;                                     \
    }                                                                          \
    [[nodiscard]] inline constexpr T##_u operator+(T a) noexcept {             \
        return static_cast<T##_u>(a);                                          \
    }                                                                          \
                                                                               \
    /* named helpers */                                                        \
    [[nodiscard]] inline constexpr bool has_flag(T mask, T flag) noexcept {    \
        return (mask & flag) == flag;                                          \
    }                                                                          \
    [[nodiscard]] inline constexpr bool has_any_of(T mask, T flag) noexcept {  \
        return (mask & flag) != static_cast<T>(0);                             \
    }





