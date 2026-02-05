#pragma once
#include <cmath>
#include <string>
#include <format>
#include <concepts>

#include "april/base/macros.hpp"
#include "april/simd/simd_traits.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/utility/debug.hpp"


namespace april::math {



    // Concepts
    template <typename T>
    concept IsScalar =  std::floating_point<T> || std::integral<T>;

    template<typename T>
    concept IsVectorSuitable = IsScalar<T> || april::simd::IsSimdType<T>;

    template <typename V>
    concept IsVectorLike = requires(V v) {
        { v.x };
        { v.y };
        { v.z };
    };



    // Forward Declarations
    template <IsVectorSuitable T, typename Scalar=T>
    struct Vec3;

    template <typename T>
    struct Vec3Proxy;

    // needed for ADL to work
    template<IsScalar T>
    AP_FORCE_INLINE T rsqrt(T val) {
        return T(1) / std::sqrt(val);
    }


    // ----------------
    // VECTOR OPS MIXIN
    // ----------------
    template <IsVectorSuitable T, typename Scalar>
    struct Vec3Ops {
        using type = T;

        // -------------------------
        // ARITHMETIC (VECTOR-VECTOR
        // -------------------------
        // vector addition
        template <IsVectorLike Other>
        Vec3<T> operator+(this const auto& self, const Other& other) noexcept {
            return {
                self.x + static_cast<T>(other.x),
                self.y + static_cast<T>(other.y),
                self.z + static_cast<T>(other.z)
            };
        }

        // vector subtraction
        template <IsVectorLike Other>
        Vec3<T> operator-(this const auto& self, const Other& other) noexcept {
            return {
                self.x - static_cast<T>(other.x),
                self.y - static_cast<T>(other.y),
                self.z - static_cast<T>(other.z)
            };
        }

        // point-wise multiplication
        template <IsVectorLike Other>
        Vec3<T> operator*(this const auto& self, const Other& other) noexcept {
            return {
                self.x * static_cast<T>(other.x),
                self.y * static_cast<T>(other.y),
                self.z * static_cast<T>(other.z)
            };
        }

        template <IsVectorLike Other>
        Vec3<T> hadamard(this const auto& self, const Other& other) noexcept {
            return self * other;
        }

        // point wise division
        template <IsVectorLike Other>
        Vec3<T> operator/(this const auto& self, const Other& other) noexcept {
            return {
                self.x / static_cast<T>(other.x),
                self.y / static_cast<T>(other.y),
                self.z / static_cast<T>(other.z)
            };
        }

        template <IsVectorLike Other>
        Vec3<T> elementwise_div(this const auto& self, const Other& other) noexcept {
            return self / other;
        }

        // unary minus
        Vec3<T> operator-(this const auto& self) noexcept {
            return {-self.x, -self.y, -self.z};
        }


        // --------------------------
        // ARITHMETIC (VECTOR-SCALAR)
        // --------------------------

        // scalar multiplication
        Vec3<T> operator*(this const auto& self, T scalar) noexcept {
            return {self.x * scalar, self.y * scalar, self.z * scalar};
        }

        // scalar division
        Vec3<T> operator/(this const auto& self, const T scalar) noexcept {
            return {self.x / scalar, self.y / scalar, self.z / scalar};
        }

        // Friend function for scalar multiplication with the scalar on the left
        template <typename S>
        requires std::convertible_to<S, T>
        friend Vec3<T> operator*(S scalar, const auto& rhs) noexcept {
            return rhs * static_cast<T>(scalar);
        }


        // --------------------
        // COMPOUND ASSIGNMENTS
        // --------------------
        template <IsVectorLike Other>
        auto& operator+=(this auto& self, const Other& rhs) noexcept {
            self.x += static_cast<T>(rhs.x);
            self.y += static_cast<T>(rhs.y);
            self.z += static_cast<T>(rhs.z);
            return self;
        }

        template <IsVectorLike Other>
        auto& operator-=(this auto& self, const Other& rhs) noexcept {
            self.x -= static_cast<T>(rhs.x);
            self.y -= static_cast<T>(rhs.y);
            self.z -= static_cast<T>(rhs.z);
            return self;
        }

        auto& operator*=(this auto& self, T scalar) noexcept {
            self.x *= scalar;
            self.y *= scalar;
            self.z *= scalar;
            return self;
        }


        // -------------------
        // GEOMETRIC FUNCTIONS
        // -------------------
        template <IsVectorLike Other>
        T dot(this const auto& self, const Other& rhs) noexcept {
            return
                self.x * static_cast<T>(rhs.x) +
                self.y * static_cast<T>(rhs.y) +
                self.z * static_cast<T>(rhs.z);
        }

        [[nodiscard]] Scalar norm_squared(this const auto& self) noexcept {
            return self.dot(self);
        }

        [[nodiscard]] Scalar norm(this const auto& self) noexcept {
            using std::sqrt;
            return sqrt(self.norm_squared());
        }

        [[nodiscard]] Scalar inv_norm(this const auto& self) noexcept {
            return rsqrt(self.norm_squared()); // compiler may optimize with fast inverse square root
        }

        [[nodiscard]] Scalar inv_norm_sq(this const auto& self) noexcept {
            return 1 / self.norm_squared(); // compiler may optimize with fast inverse square root
        }



        // -------------------
        // ORDERING & EQUALITY
        // -------------------
        // v < u iff for all v_i: v_i < u_i
        // in other words v smaller than u if every element in v is smaller than the corresponding element in u

        template <IsVectorLike Other>
        bool operator==(this const auto& self, const Other& other) noexcept {
            return
                self.x == static_cast<T>(other.x) &&
                self.y == static_cast<T>(other.y) &&
                self.z == static_cast<T>(other.z);
        }

        template <IsVectorLike Other>
        bool operator<=(this const auto& self, const Other& other) noexcept {
            return self.x <= static_cast<T>(other.x) &&
                   self.y <= static_cast<T>(other.y) &&
                   self.z <= static_cast<T>(other.z);
        }

        template <IsVectorLike Other>
        bool operator>=(this const auto& self, const Other& other) noexcept {
            return self.x >= static_cast<T>(other.x) &&
                   self.y >= static_cast<T>(other.y) &&
                   self.z >= static_cast<T>(other.z);
        }

        template <IsVectorLike Other>
        bool operator<(this const auto& self, const Other& other) noexcept {
            return self.x < static_cast<T>(other.x) &&
                   self.y < static_cast<T>(other.y) &&
                   self.z < static_cast<T>(other.z);
        }

        template <IsVectorLike Other>
        bool operator>(this const auto& self, const Other& other) noexcept {
            return self.x > static_cast<T>(other.x) &&
                   self.y > static_cast<T>(other.y) &&
                   self.z > static_cast<T>(other.z);
        }


        // ---------
        // ACCESSORS
        // ---------
        // Access component by index: 0 for x, 1 for y, 2 for z.
        // decltype(auto) preserves references
        decltype(auto) operator[](this auto&& self, const int index) noexcept {
            AP_ASSERT(index >= 0 && index < 3, "Index out of bounds");
            if (index == 0) return (self.x); // Parentheses matter for decltype(auto) on members
            if (index == 1) return (self.y);
            return (self.z);
        }

        T max(this const auto& self) noexcept {
            return std::max(self.x, std::max(self.y, self.z));
        }

        T min(this const auto& self) noexcept {
            return std::min(self.x, std::min(self.y, self.z));
        }


        //-----------------
        // LOGIC PREDICATES
        // ----------------
        template <typename Predicate>
        bool any(this const auto& self, Predicate predicate) {
            return predicate(self.x) || predicate(self.y) || predicate(self.z);
        }

        template <typename Predicate>
        bool all(this const auto& self, Predicate predicate) {
            return predicate(self.x) && predicate(self.y) && predicate(self.z);
        }

        // debug print
        [[nodiscard]] std::string to_string(this const auto& self) {
            return std::format("{{{}, {}, {}}}", self.x, self.y, self.z);
        }
    };




    // ---------------
    // VEC3 DEFINITION
    // ---------------
    template <IsVectorSuitable T, typename Scalar>
    struct Vec3 : Vec3Ops<T, Scalar> {
        T x, y, z;

        Vec3() : x(0), y(0), z(0) {}
        Vec3(T x, T y, T z) : x(x), y(y), z(z) {}
        explicit Vec3(T v) : x(v), y(v), z(v) {}

        template <IsVectorLike Other>
        Vec3(const Other& p) : x(static_cast<T>(p.x)), y(static_cast<T>(p.y)), z(static_cast<T>(p.z)) {}
    };



    // ------------
    // VEC3 POINTER
    // ------------
    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3Ptr {
        T * AP_RESTRICT x = nullptr;
        T * AP_RESTRICT y = nullptr;
        T * AP_RESTRICT z = nullptr;

        Vec3Ptr() = default;

        template <typename U>
        requires std::convertible_to<U, T>
        Vec3Ptr(Vec3<U>* other)
           : x(&other->x), y(&other->y), z(&other->z) {
            AP_ASSERT(x != y && y != z &&  z!= x, "x y z pointers do not point to different addresses");
        }

        template <typename U>
        requires std::convertible_to<const U*, T*>
        Vec3Ptr(const Vec3<U>* AP_RESTRICT other)
           : x(&other->x), y(&other->y), z(&other->z) {
            AP_ASSERT(x != y && y != z &&  z!= x, "x y z pointers do not point to different addresses");
        }

        Vec3Ptr(T& AP_RESTRICT x_ref, T& AP_RESTRICT y_ref, T& AP_RESTRICT z_ref)
           : x(&x_ref), y(&y_ref), z(&z_ref) {
            AP_ASSERT(x != y && y != z &&  z!= x, "x y z pointers do not point to different addresses");
        }

        Vec3Ptr(T* AP_RESTRICT x_ptr, T* AP_RESTRICT y_ptr, T* AP_RESTRICT z_ptr)
            : x(x_ptr), y(y_ptr), z(z_ptr) {
            AP_ASSERT(x != y && y != z &&  z!= x, "x y z pointers do not point to different addresses");
        }

        template <typename U>
        requires std::convertible_to<U*, T*>
        Vec3Ptr(const Vec3Ptr<U>& other)
           : x(other.x), y(other.y), z(other.z) {
            AP_ASSERT(x != y && y != z &&  z!= x, "x y z pointers do not point to different addresses");
        }

        template <typename U>
        requires std::convertible_to<U*, T*>
        Vec3Ptr& operator=(const Vec3Ptr<U>& other) {
            x = other.x;
            y = other.y;
            z = other.z;
            return *this;
        }

        Vec3Proxy<T> operator*() const {
            return Vec3Proxy<T>(*x, *y, *z);
        }
    };



    // ----------
    // VEC3 PROXY
    // ----------
    template <typename T>
    struct Vec3Proxy : Vec3Ops<T, double> {
        T& AP_RESTRICT x;
        T& AP_RESTRICT y;
        T& AP_RESTRICT z;

        Vec3Proxy(const Vec3Proxy&) = default;

        Vec3Proxy(T& x_ref, T& y_ref, T& z_ref)
            : x(x_ref), y(y_ref), z(z_ref) {}

        template<typename U>
        requires std::convertible_to<U, T>
        explicit Vec3Proxy(Vec3<U> & other)
            : x(other.x), y(other.y), z(other.z) {}

        template <typename U>
        requires std::convertible_to<U&, T&>
        explicit Vec3Proxy(const Vec3Proxy<U>& other)
            : x(other.x), y(other.y), z(other.z) {}

        Vec3Proxy& operator=(const Vec3<T>& rhs) {
            x = rhs.x; y = rhs.y; z = rhs.z;
            return *this;
        }

        Vec3Proxy& operator=(const Vec3Proxy& rhs) {
            x = rhs.x; y = rhs.y; z = rhs.z;
            return *this;
        }

        // implicit conversion to Value
        operator Vec3<T>() const { return Vec3<T>(x, y, z); }
    };

    // specialization for packed types
    template <simd::IsSimdType T>
    struct Vec3Proxy<T> : Vec3Ops<T, T> {
        using Ref = simd::PackedRef<typename T::value_type, T>;
        using Scalar = T::value_type;

        // Members are "Reference Wrappers", not C++ references
        Ref x;
        Ref y;
        Ref z;
        Vec3Proxy(const Vec3Proxy&) = default;

        template<typename U>
        Vec3Proxy(const Vec3Ptr<U> & ptr): x(ptr.x), y(ptr.y), z(ptr.z) {}


        Vec3Proxy(T& x_ref, T& y_ref, T& z_ref)
            : x(x_ref), y(y_ref), z(z_ref) {}

        template<typename U>
        requires std::convertible_to<U, T>
        explicit Vec3Proxy(Vec3<U> & other)
            : x(other.x), y(other.y), z(other.z) {}

        template <typename U>
        requires std::convertible_to<U&, T&>
        explicit Vec3Proxy(const Vec3Proxy<U>& other)
            : x(other.x), y(other.y), z(other.z) {}

        // Constructor from POINTERS (SoA style)
        // This is crucial for SIMD iterators!
        Vec3Proxy(Scalar* ptr_x, Scalar* ptr_y, Scalar* ptr_z)
            : x(ptr_x), y(ptr_y), z(ptr_z) {}

        // Assignment from Value (Vec3<Packed>)
        // "p.pos = result_vec;"
        Vec3Proxy& operator=(const Vec3<T>& rhs) {
            x = rhs.x; // Calls PackedRef::operator=(Packed) -> Stores to memory
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        // Assignment from other Proxy (Copy Memory to Memory)
        // "p.pos = p.old_pos;"
        Vec3Proxy& operator=(const Vec3Proxy& rhs) {
            if (this != &rhs) {
                x = rhs.x; // Calls PackedRef::operator=(PackedRef)
                y = rhs.y;
                z = rhs.z;
            }
            return *this;
        }

        // Implicit conversion to Value
        operator Vec3<T>() const { return Vec3<T>(x, y, z); }
    };



} // namespace april::utils