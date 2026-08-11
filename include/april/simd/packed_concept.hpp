#pragma once
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <utility>


namespace april::simd {

    template<typename T, bool = std::is_enum_v<std::remove_cv_t<T>>>
    struct packed_storage {
        using type = std::remove_cv_t<T>;
    };

    template<typename T>
    struct packed_storage<T, true> {
        using type = std::underlying_type_t<std::remove_cv_t<T>>;
    };

    template<typename T>
    using packed_storage_t = packed_storage<T>::type;



    template<size_t N>
    struct ByteOffsets {
        using value_type = std::ptrdiff_t;

        std::array<value_type, N> values;

        template<std::integral... T>
            requires (sizeof...(T) == N)
        constexpr ByteOffsets(T... offsets) noexcept
            : values{static_cast<value_type>(offsets)...}
        {}

        template<std::integral T>
        constexpr ByteOffsets(const std::array<T, N>& offsets) noexcept {
            for (size_t i = 0; i < N; ++i)
                values[i] = static_cast<value_type>(offsets[i]);
        }
    };

    template<typename... T>
    ByteOffsets(T...) -> ByteOffsets<sizeof...(T)>;

    template<std::integral T, size_t N>
    ByteOffsets(const std::array<T, N>&) -> ByteOffsets<N>;


    template<typename T>
    concept IsPackableValue =
        std::is_arithmetic_v<std::remove_cv_t<T>> ||
        std::is_enum_v<std::remove_cv_t<T>>;


    // mask concept
    template<typename M>
    concept HasMaskReductions = requires(const std::remove_cvref_t<M>& mask) {
        { std::remove_cvref_t<M>::all(mask) }  -> std::same_as<bool>;
        { std::remove_cvref_t<M>::any(mask) }  -> std::same_as<bool>;
        { std::remove_cvref_t<M>::none(mask) } -> std::same_as<bool>;
    };

    template<typename T>
    concept IsSimdMaskImpl = HasMaskReductions<T> &&
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


    template<typename T>
    concept HasCommonMathFunctions = requires(
        const T a, const T b, const T c,
        const typename T::mask_type mask
    ) {
        // Selection
        { T::select(mask, a, b) } -> std::same_as<T>;

        // Basic
        { T::abs(a) }       -> std::same_as<T>;
        { T::min(a, b) }    -> std::same_as<T>;
        { T::max(a, b) }    -> std::same_as<T>;
        { T::clamp(a, b, c) } -> std::same_as<T>;

        // Rounding: identity operations for integer packs
        { T::floor(a) }     -> std::same_as<T>;
        { T::ceil(a) }      -> std::same_as<T>;
        { T::round(a) }     -> std::same_as<T>;
        { T::trunc(a) }     -> std::same_as<T>;
        { T::nearbyint(a) } -> std::same_as<T>;

        // Numeric
        { T::fma(a, b, c) } -> std::same_as<T>;

        // Classification
        { T::isnan(a) }    -> std::same_as<typename T::mask_type>;
        { T::isinf(a) }    -> std::same_as<typename T::mask_type>;
        { T::isfinite(a) } -> std::same_as<typename T::mask_type>;
        { T::signbit(a) }  -> std::same_as<typename T::mask_type>;
    };

    template<typename T>
    concept HasFloatingMathFunctions =
    std::floating_point<typename T::storage_type> &&
    requires(const T a, const T b) {
        // Roots and powers
        { T::sqrt(a) }     -> std::same_as<T>;
        { T::rsqrt(a) }    -> std::same_as<T>;
        { T::cbrt(a) }     -> std::same_as<T>;
        { T::hypot(a, b) } -> std::same_as<T>;
        { T::pow(a, b) }   -> std::same_as<T>;

        // Exponential and logarithmic
        { T::exp(a) }   -> std::same_as<T>;
        { T::exp2(a) }  -> std::same_as<T>;
        { T::expm1(a) } -> std::same_as<T>;
        { T::log(a) }   -> std::same_as<T>;
        { T::ln(a) }    -> std::same_as<T>;
        { T::log2(a) }  -> std::same_as<T>;
        { T::log10(a) } -> std::same_as<T>;
        { T::log1p(a) } -> std::same_as<T>;

        // Trigonometric
        { T::sin(a) }       -> std::same_as<T>;
        { T::cos(a) }       -> std::same_as<T>;
        { T::sincos(a) }    -> std::same_as<std::pair<T, T>>;
        { T::tan(a) }       -> std::same_as<T>;
        { T::asin(a) }      -> std::same_as<T>;
        { T::acos(a) }      -> std::same_as<T>;
        { T::atan(a) }      -> std::same_as<T>;
        { T::atan2(a, b) }  -> std::same_as<T>;

        // Hyperbolic
        { T::sinh(a) }  -> std::same_as<T>;
        { T::cosh(a) }  -> std::same_as<T>;
        { T::tanh(a) }  -> std::same_as<T>;
        { T::asinh(a) } -> std::same_as<T>;
        { T::acosh(a) } -> std::same_as<T>;
        { T::atanh(a) } -> std::same_as<T>;

        // Floating-point numeric operations
        { T::fmod(a, b) }      -> std::same_as<T>;
        { T::remainder(a, b) } -> std::same_as<T>;
        { T::copysign(a, b) }  -> std::same_as<T>;
    };

    template<typename T>
    concept HasMathFunctions =
        HasCommonMathFunctions<T> &&
        (!std::floating_point<typename T::storage_type> || HasFloatingMathFunctions<T>);


    template<typename T>
    concept HasReductionsOps = requires(const T ct) {
        { ct.reduce_add() } -> std::same_as<typename T::storage_type>;
        { ct.reduce_min() } -> std::same_as<typename T::storage_type>;
        { ct.reduce_max() } -> std::same_as<typename T::storage_type>;
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
        && HasScalarMixedOps<T, typename T::storage_type>;

    // The Main Concept
    template<typename T>
    concept IsSimdTypeImpl = requires(
        T t,
        const T ct,
        typename T::value_type scalar,
        const typename T::value_type* ptr,
        ByteOffsets<T::size()> byte_offsets,
        std::array<const typename T::value_type*, T::size()> gather_ptrs,
        std::array<typename T::value_type*, T::size()> scatter_ptrs,
        std::ptrdiff_t byte_stride
    ) {
        typename T::value_type;
        typename T::mask_type;
        typename T::storage_type;

        requires std::is_arithmetic_v<typename T::storage_type>;
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
        // base + lane * byte_stride
        { T::gather_strided(ptr, byte_stride) }-> std::same_as<T>;
        // Compile-time byte stride
        { T::template gather_strided<sizeof(typename T::value_type)>(ptr) }-> std::same_as<T>;
        // Arbitrary byte offsets: base + byte_offsets[lane]
        { T::gather(ptr, byte_offsets) }-> std::same_as<T>;
        // Completely arbitrary addresses
        { T::gather(gather_ptrs) }-> std::same_as<T>;

        // storing
        { ct.store(const_cast<typename T::value_type*>(ptr)) } -> std::same_as<void>;
        { ct.store_aligned(const_cast<typename T::value_type*>(ptr)) } -> std::same_as<void>;
        { ct.store_unaligned(const_cast<typename T::value_type*>(ptr)) } -> std::same_as<void>;

        // scattering
        // Regular byte stride
        { ct.scatter_strided(const_cast<T::value_type*>(ptr),byte_stride) } -> std::same_as<void>;
        // Compile-time byte stride
        { ct.template scatter_strided<sizeof(typename T::value_type)>(const_cast<T::value_type*>(ptr)) } -> std::same_as<void>;
        // Arbitrary byte offsets
        { ct.scatter(const_cast<T::value_type*>(ptr), byte_offsets) } -> std::same_as<void>;
        // Completely arbitrary addresses
        { ct.scatter(scatter_ptrs) } -> std::same_as<void>;


        // masking
        { t == t } -> std::same_as<typename T::mask_type>;
        { T::select(t == t, t, t) } -> std::same_as<T>;

        // permutations
        { ct.rotate_left() } -> std::same_as<T>;
        { ct.rotate_right() } -> std::same_as<T>;
        { ct.template rotate_left<2>() } -> std::same_as<T>;
        { ct.template rotate_right<2>() } -> std::same_as<T>;
        { ct.template permute<0>() } -> std::same_as<T>;

        // bitwise requirement (only enforced if the scalar type is an integer)
        requires (std::is_integral_v<typename T::storage_type> ? HasBitwiseOps<T> : true);

    } && HasSimdOps<T>;

    template<typename T>
    concept IsSimdType = IsSimdTypeImpl<std::remove_cvref_t<T>>;
}


