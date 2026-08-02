#pragma once
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>


namespace april::simd {

    // mask concept
    template<typename T>
    concept IsSimdMaskImpl =
        std::default_initializable<T> &&
        std::copyable<T> &&
        requires(T m, T m2, const T cm, bool* ptr, const bool* cptr) {
        // Type construction and width
        { T{false} } -> std::same_as<T>;
        { T::size() } -> std::same_as<std::size_t>;

        // Memory
        { T::load(cptr) }           -> std::same_as<T>;
        { T::load_aligned(cptr) }   -> std::same_as<T>;
        { T::load_unaligned(cptr) } -> std::same_as<T>;

        { cm.store(ptr) }           -> std::same_as<void>;
        { cm.store_aligned(ptr) }   -> std::same_as<void>;
        { cm.store_unaligned(ptr) } -> std::same_as<void>;

        // Bitmask import/export
        { cm.to_bitmask() } -> std::same_as<std::uint64_t>;
        { T::from_bitmask(std::uint64_t{}) } -> std::same_as<T>;
        { cm.to_array() } -> std::same_as<std::array<bool, T::size()>>;
        { cm.to_string() } -> std::same_as<std::string>;

        // Reductions
        { all(cm) }  -> std::same_as<bool>;
        { any(cm) }  -> std::same_as<bool>;
        { none(cm) } -> std::same_as<bool>;

        // Logical operators
        { !cm }       -> std::same_as<T>;
        { cm && m2 }  -> std::same_as<T>;
        { cm || m2 }  -> std::same_as<T>;

        // Bitwise operators
        { cm & m2 } -> std::same_as<T>;
        { cm | m2 } -> std::same_as<T>;
        { cm ^ m2 } -> std::same_as<T>;
        { ~cm }     -> std::same_as<T>;

        // Lane-wise comparisons
        { cm == m2 } -> std::same_as<T>;
        { cm != m2 } -> std::same_as<T>;

        // Mutating lane rotations
        { m.rotate_left() }             -> std::same_as<void>;
        { m.rotate_right() }            -> std::same_as<void>;
        { m.template rotate_left<2>() } -> std::same_as<void>;
        { m.template rotate_right<2>() }-> std::same_as<void>;
        };

    template<typename T>
    concept IsSimdMask = IsSimdMaskImpl<std::remove_cvref_t<T>>;

    // check if all usual arithmetic ops exist
    template<typename T>
    concept HasArithmeticOps = requires(T a, T b) {
        { + a } -> std::same_as<T>;
        { - a } -> std::same_as<T>;
        { a + b } -> std::same_as<T>;
        { a - b } -> std::same_as<T>;
        { a * b } -> std::same_as<T>;
        { a / b } -> std::same_as<T>;
        { a += b } -> std::same_as<T&>;
        { a -= b } -> std::same_as<T&>;
        { a *= b } -> std::same_as<T&>;
        { a /= b } -> std::same_as<T&>;
    };

    // check if comparator ops exist
    template<typename T>
    concept HasComparisonOps = requires(T a, T b) {
        typename T::mask_type;
        requires IsSimdMask<typename T::mask_type>;

        { a == b } -> std::same_as<typename T::mask_type>;
        { a != b } -> std::same_as<typename T::mask_type>;
        { a < b }  -> std::same_as<typename T::mask_type>;
        { a <= b } -> std::same_as<typename T::mask_type>;
        { a > b }  -> std::same_as<typename T::mask_type>;
        { a >= b } -> std::same_as<typename T::mask_type>;
    };

    template<typename T, typename Scalar>
    concept HasScalarMixedOps = requires(T t, Scalar s) {
        // Arithmetic (Left & Right)
        { t + s } -> std::same_as<T>;
        { s + t } -> std::same_as<T>;
        { t - s } -> std::same_as<T>;
        { s - t } -> std::same_as<T>;
        { t * s } -> std::same_as<T>;
        { s * t } -> std::same_as<T>;
        { t / s } -> std::same_as<T>;
        { s / t } -> std::same_as<T>;

        // Compound (Vector on LHS only)
        { t += s } -> std::same_as<T&>;
        { t -= s } -> std::same_as<T&>;
        { t *= s } -> std::same_as<T&>;
        { t /= s } -> std::same_as<T&>;

        // Comparison (Left & Right)
        { t == s } -> std::same_as<typename T::mask_type>;
        { s == t } -> std::same_as<typename T::mask_type>;
        { t != s } -> std::same_as<typename T::mask_type>;
        { s != t } -> std::same_as<typename T::mask_type>;
        { t < s }  -> std::same_as<typename T::mask_type>;
        { s < t }  -> std::same_as<typename T::mask_type>;
        { t <= s } -> std::same_as<typename T::mask_type>;
        { s <= t } -> std::same_as<typename T::mask_type>;
        { t > s }  -> std::same_as<typename T::mask_type>;
        { s > t }  -> std::same_as<typename T::mask_type>;
        { t >= s } -> std::same_as<typename T::mask_type>;
        { s >= t } -> std::same_as<typename T::mask_type>;
    };


    // check for free functions
    template<typename T>
    concept HasMathFunctions = requires(T a, T b, T c) {
        { sqrt(a) }     -> std::same_as<T>;
        { rsqrt(a) }    -> std::same_as<T>;
        { abs(a) }      -> std::same_as<T>;
        { min(a, b) }   -> std::same_as<T>;
        { max(a, b) }   -> std::same_as<T>;
        { fma(a, b, c) }-> std::same_as<T>;
        { round(a) }    -> std::same_as<T>;
        { floor(a) }    -> std::same_as<T>;
        { ceil(a) }     -> std::same_as<T>;
    };

    template<typename T>
    concept HasReductionsOps = requires(const T ct) {
        { ct.reduce_add() } -> std::same_as<typename T::value_type>;
        { ct.reduce_min() } -> std::same_as<typename T::value_type>;
        { ct.reduce_max() } -> std::same_as<typename T::value_type>;
    };

    template<typename T> // only for integers
    concept HasBitwiseOps = requires(T a, T b) {
        { ~a } -> std::same_as<T>;
        { a & b } -> std::same_as<T>;
        { a | b } -> std::same_as<T>;
        { a ^ b } -> std::same_as<T>;
        { a &= b } -> std::same_as<T&>;
        { a |= b } -> std::same_as<T&>;
        { a ^= b } -> std::same_as<T&>;
    };

    template<typename T>
    concept HasSimdOps =
        HasArithmeticOps<T>
        && HasComparisonOps<T>
        && HasMathFunctions<T>
        && HasReductionsOps<T>
        && HasScalarMixedOps<T, float>
        && HasScalarMixedOps<T, double>
        && HasScalarMixedOps<T, long double>;

    // The Main Concept
    template<typename T>
    concept IsSimdTypeImpl = requires(T t, const T ct, typename T::value_type scalar, const typename T::value_type* ptr) {
        typename T::value_type;
        typename T::mask_type;

        requires std::is_arithmetic_v<typename T::value_type>;
        requires IsSimdMask<typename T::mask_type>;

        { T::size() } -> std::convertible_to<std::size_t>;
        { ct.to_string() } -> std::convertible_to<std::string>;

        T();
        T(scalar); // Broadcast
        { t = scalar } -> std::same_as<T&>; // Scalar Assignment (Broadcast)

        // loading
        { T::load(ptr) } -> std::same_as<T>;
        { T::load_aligned(ptr) } -> std::same_as<T>;
        { T::load_unaligned(ptr) } -> std::same_as<T>;

        // gathering
        // Indirect load via pointer array
        { T::gather(static_cast<const T::value_type* const*>(nullptr)) } -> std::same_as<T>;
        // Indirect load via offsets (using a vector of the same width for indices)
        { T::gather(ptr, t) } -> std::same_as<T>;

        // storing
        { ct.store(const_cast<T::value_type*>(ptr)) } -> std::same_as<void>;
        { ct.store_aligned(const_cast<T::value_type*>(ptr)) } -> std::same_as<void>;
        { ct.store_unaligned(const_cast<T::value_type*>(ptr)) } -> std::same_as<void>;

        // scattering
        // Indirect store via offsets
        { ct.scatter(const_cast<T::value_type*>(ptr), t) } -> std::same_as<void>;

        // masking
        { t == t } -> std::same_as<typename T::mask_type>;
        { select(t == t, t, t) } -> std::same_as<T>;

        // permutations
        { ct.rotate_left() } -> std::same_as<T>;
        { ct.rotate_right() } -> std::same_as<T>;
        { ct.template rotate_left<2>() } -> std::same_as<T>;
        { ct.template rotate_right<2>() } -> std::same_as<T>;
        { ct.template permute<0>() } -> std::same_as<T>;

        // bitwise requirement (only enforced if the scalar type is an integer)
        requires (std::is_integral_v<typename T::value_type> ? HasBitwiseOps<T> : true);

    } && HasSimdOps<T>;

    template<typename T>
    concept IsSimdType = IsSimdTypeImpl<std::remove_cvref_t<T>>;
}


