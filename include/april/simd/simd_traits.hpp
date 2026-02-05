#pragma once
#include <concepts>
#include <cstddef>
#include <string>


namespace april::simd {

    // mask concept
    template<typename T>
    concept IsSimdMask = requires(T m, T m2) {
        // Reductions
        { all(m) } -> std::same_as<bool>;
        { any(m) } -> std::same_as<bool>;

        // Logical Operators
        { !m }       -> std::same_as<T>;
        { m && m2 }  -> std::same_as<T>;
        { m || m2 }  -> std::same_as<T>;

        // Comparisons (Mask == Mask)
        { m == m2 } -> std::same_as<T>;
        { m != m2 } -> std::same_as<T>;
    };

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
        { a == b };
        { a != b };
        { a < b };
        { a <= b };
        { a > b };
        { a >= b };
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
        { t == s }; { s == t };
        { t != s }; { s != t };
        { t < s };  { s < t };
        { t <= s }; { s <= t };
        { t > s };  { s > t };
        { t >= s }; { s >= t };
    };


    // check for free functions
    template<typename T>
    concept HasMathFunctions = requires(T a, T b, T c) {
        { sqrt(a) } -> std::same_as<T>;
        { rsqrt(a) } -> std::same_as<T>;
        { abs(a) } -> std::same_as<T>;
        { min(a, b) } -> std::same_as<T>;
        { max(a, b) } -> std::same_as<T>;
        { fma(a, b, c) } -> std::same_as<T>;
    };


    template<typename T>
    concept HasSimdOps =
        HasArithmeticOps<T>
        && HasComparisonOps<T>
        && HasMathFunctions<T>
        && HasScalarMixedOps<T, float>
        && HasScalarMixedOps<T, double>
        && HasScalarMixedOps<T, long double>;


    // The Main Concept
    template<typename T>
    concept IsSimdTypeImpl = requires(T t, const T ct, typename T::value_type scalar, const typename T::value_type* ptr) {
        typename T::value_type;
        requires (std::floating_point<typename T::value_type>);
        { T::size() } -> std::convertible_to<std::size_t>;
        { ct.to_string() } -> std::convertible_to<std::string>;

        T();
        T(scalar); // Broadcast

        // loading
        { T::load(ptr) } -> std::same_as<T>;
        { T::load_aligned(ptr) } -> std::same_as<T>;
        { T::load_unaligned(ptr) } -> std::same_as<T>;
        { T::gather(static_cast<const T::value_type* const*>(nullptr)) } -> std::same_as<T>;

        // storing
        ct.store(const_cast<T::value_type*>(ptr));
        ct.store_aligned(const_cast<T::value_type*>(ptr));
        ct.store_unaligned(const_cast<T::value_type*>(ptr));

        // permutations
        { ct.rotate_left() } -> std::same_as<T>;
        { ct.rotate_right() } -> std::same_as<T>;
        { ct.template rotate_left<2>() } -> std::same_as<T>;
        { ct.template rotate_right<2>() } -> std::same_as<T>;
        { ct.template permute<0>() } -> std::same_as<T>;
    } && HasSimdOps<T>;

    template<typename T>
    concept IsSimdType = IsSimdTypeImpl<std::remove_cvref_t<T>>;
}