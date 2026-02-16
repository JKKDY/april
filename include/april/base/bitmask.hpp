#pragma once

#include <type_traits>

namespace april {

    namespace bitmask::internal {
        // Opt-in Toggle
        template<typename T>
        inline constexpr bool enable_bitmask_operators = false;

        // BitMask concept
        template<typename T>
        concept BitmaskEnum = std::is_enum_v<T> && enable_bitmask_operators<T>;
    }

    // return true if flag is a subset mask
    template<bitmask::internal::BitmaskEnum T>
    [[nodiscard]] constexpr bool has_flag(T mask, T flag) noexcept {
        using U = std::underlying_type_t<T>;
        return (static_cast<U>(mask) & static_cast<U>(flag)) == static_cast<U>(flag);
    }

    // return true if flag and mask have a non-zero intersection
    template<bitmask::internal::BitmaskEnum T>
    [[nodiscard]] constexpr bool has_any_of(T mask, T flag) noexcept {
        using U = std::underlying_type_t<T>;
        return (static_cast<U>(mask) & static_cast<U>(flag)) != 0;
    }
}

#define AP_ENABLE_BITMASK_OPERATORS(T) \
    /* Bitwise OR */ \
    inline constexpr T operator|(T lhs, T rhs) noexcept { \
        return static_cast<T>(static_cast<std::underlying_type_t<T>>(lhs) | \
                              static_cast<std::underlying_type_t<T>>(rhs)); \
    } \
    inline constexpr T& operator|=(T& lhs, T rhs) noexcept { \
        return lhs = lhs | rhs; \
    } \
    /* Bitwise AND */ \
    inline constexpr T operator&(T lhs, T rhs) noexcept { \
        return static_cast<T>(static_cast<std::underlying_type_t<T>>(lhs) & \
                              static_cast<std::underlying_type_t<T>>(rhs)); \
    } \
    inline constexpr T& operator&=(T& lhs, T rhs) noexcept { \
        return lhs = lhs & rhs; \
    } \
    /* Bitwise XOR */ \
    inline constexpr T operator^(T lhs, T rhs) noexcept { \
        return static_cast<T>(static_cast<std::underlying_type_t<T>>(lhs) ^ \
                              static_cast<std::underlying_type_t<T>>(rhs)); \
    } \
    inline constexpr T& operator^=(T& lhs, T rhs) noexcept { \
        return lhs = lhs ^ rhs; \
    } \
    /* Bitwise NOT */ \
    inline constexpr T operator~(T val) noexcept { \
        return static_cast<T>(~static_cast<std::underlying_type_t<T>>(val)); \
    } \
    /* Set Difference (Remove flag) */ \
    inline constexpr T operator-(T lhs, T rhs) noexcept { \
        return lhs & ~rhs; \
    } \
    inline constexpr T& operator-=(T& lhs, T rhs) noexcept { \
        return lhs = lhs & ~rhs; \
    } \
    /* Utility: Logical check (if (!flag)) */ \
    inline constexpr bool operator!(T val) noexcept { \
        return static_cast<std::underlying_type_t<T>>(val) == 0; \
    } \
    /* Utility: Unary plus for underlying cast */ \
    inline constexpr std::underlying_type_t<T> operator+(T val) noexcept { \
        return static_cast<std::underlying_type_t<T>>(val); \
    }






