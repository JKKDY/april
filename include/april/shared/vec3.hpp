#pragma once
#include <cmath>
#include <string>
#include <format>
#include <concepts>
#include <functional>

#include "april/shared/debug.hpp"


namespace april::utils {


    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3;

    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3Proxy;


    template <typename V>
    concept IsVectorLike = requires(V v) {
        { v.x };
        { v.y };
        { v.z };
    };

    template <typename S, typename T>
    concept IsScalar = std::convertible_to<S, T>; // Or std::arithmetic<S>


    template <typename T>
    requires (std::integral<std::remove_cvref_t<T>> || std::floating_point<std::remove_cvref_t<T>>)
    struct Vec3Ops {
        using type = T;

        // --------------
        // ARITHMETIC OPS
        // --------------
        // vector addition
        template <IsVectorLike Other>
        Vec3<T> operator+(this const auto& self, const Other& other) noexcept {
            return {self.x + other.x, self.y + other.y, self.z + other.z};
        }

        // vector subtraction
        template <IsVectorLike Other>
        Vec3<T> operator-(this const auto& self, const Other& other) noexcept {
            return {self.x - other.x, self.y - other.y, self.z - other.z};
        }

        // point-wise multiplication
        template <IsVectorLike Other>
        Vec3<T> operator*(this const auto& self, const Other& other) noexcept {
            return {self.x * other.x, self.y * other.y, self.z * other.z};
        }

        template <IsVectorLike Other>
        Vec3<T> hadamard(this const auto& self, const Other& other) noexcept {
            return self * other;
        }

        // point wise division
        template <IsVectorLike Other>
        Vec3<T> operator/(this const auto& self, const Other& other) noexcept {
            return {self.x / other.x, self.y / other.y, self.z / other.z};
        }

        template <IsVectorLike Other>
        Vec3<T> elementwise_div(this const auto& self, const Other& other) noexcept {
            return self / other;
        }

        // unary minus
        Vec3<T> operator-(this const auto& self) noexcept {
            return {-self.x, -self.y, -self.z};
        }

        // scalar multiplication
        Vec3<T> operator*(this const auto& self, T scalar) noexcept {
            return {self.x * scalar, self.y * scalar, self.z * scalar};
        }

        // scalar division
        Vec3<T> operator/(this const auto& self, const T scalar) noexcept {
            return {self.x / scalar, self.y / scalar, self.z / scalar};
        }

        // Friend function for scalar multiplication with the scalar on the left
        template <typename Scalar>
        requires std::convertible_to<Scalar, T>
        friend Vec3<T> operator*(Scalar scalar, const auto& rhs) noexcept {
            return rhs * static_cast<T>(scalar);
        }


        // --------------------
        // COMPOUND ASSIGNMENTS
        // --------------------
        template <IsVectorLike Other>
        auto& operator+=(this auto& self, const Other& rhs) noexcept {
            self.x += rhs.x;
            self.y += rhs.y;
            self.z += rhs.z;
            return self;
        }

        template <IsVectorLike Other>
        auto& operator-=(this auto& self, const Other& rhs) noexcept {
            self.x -= rhs.x;
            self.y -= rhs.y;
            self.z -= rhs.z;
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
            return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z;
        }

        [[nodiscard]] double norm_squared(this const auto& self) noexcept {
            return self.x * self.x + self.y * self.y + self.z * self.z;
        }

        [[nodiscard]] double norm(this const auto& self) noexcept {
            return std::sqrt(self.norm_squared());
        }

        [[nodiscard]] double inv_norm(this const auto& self) noexcept {
            return 1 / std::sqrt(self.norm_squared()); // compiler may optimize with fast inverse square root
        }


        // -------------------
        // ORDERING & EQUALITY
        // -------------------
        // v < u iff for all v_i: v_i < u_i
        // in other words v smaller than u if every element in v is smaller than the corresponding element in u

        template <IsVectorLike Other>
        bool operator==(this const auto& self, const Other& other) noexcept {
            return self.x == other.x && self.y == other.y && self.z == other.z;
        }

        template <IsVectorLike Other>
        bool operator<=(this const auto& self, const Other& other) noexcept {
            return self.x <= other.x && self.y <= other.y && self.z <= other.z;
        }

        template <IsVectorLike Other>
        bool operator>=(this const auto& self, const Other& other) noexcept {
            return self.x >= other.x && self.y >= other.y && self.z >= other.z;
        }

        template <IsVectorLike Other>
        bool operator<(this const auto& self, const Other& other) noexcept {
            return self.x < other.x && self.y < other.y && self.z < other.z;
        }

        template <IsVectorLike Other>
        bool operator>(this const auto& self, const Other& other) noexcept {
            return self.x > other.x && self.y > other.y && self.z > other.z;
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


        [[nodiscard]] std::string to_string(this const auto& self) {
            return std::format("{{{}, {}, {}}}", self.x, self.y, self.z);
        }
    };



    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3 : Vec3Ops<T> {
        T x, y, z;

        Vec3() : x(0), y(0), z(0) {}
        Vec3(T x, T y, T z) : x(x), y(y), z(z) {}
        explicit Vec3(T v) : x(v), y(v), z(v) {}

        template <IsVectorLike Other>
        Vec3(const Other& p) : x(p.x), y(p.y), z(p.z) {}
    };



    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3Ptr {
        T * x = nullptr;
        T * y = nullptr;
        T * z = nullptr;

        Vec3Ptr() = default;

        template <typename U>
        requires std::convertible_to<U, T>
        Vec3Ptr(Vec3<U>* other)
           : x(&other->x), y(&other->y), z(&other->z) {}

        template <typename U>
        requires std::convertible_to<const U*, T*>
        Vec3Ptr(const Vec3<U>* other)
           : x(&other->x), y(&other->y), z(&other->z) {}

        Vec3Ptr(T& x_ref, T& y_ref, T& z_ref)
           : x(&x_ref), y(&y_ref), z(&z_ref) {}

        Vec3Ptr(T* x_ptr, T* y_ptr, T* z_ptr)
            : x(x_ptr), y(y_ptr), z(z_ptr) {}

        template <typename U>
        requires std::convertible_to<U*, T*>
        Vec3Ptr(const Vec3Ptr<U>& other)
           : x(other.x), y(other.y), z(other.z) {}

        template <typename U>
        requires std::convertible_to<U*, T*>
        Vec3Ptr& operator=(const Vec3Ptr<U>& other) {
            x = other.x;
            y = other.y;
            z = other.z;
            return *this;
        }

        // Vec3Ptr(T* x_ptr, T* y_ptr, T* z_ptr)
        //   : x(x_ptr), y(y_ptr), z(z_ptr) {}

        Vec3Proxy<T> operator*() const {
            return Vec3Proxy<T>(*x, *y, *z);
        }
    };


    template <typename T> requires std::integral<T> || std::floating_point<T>
    struct Vec3Proxy : Vec3Ops<T> {
        T& x;
        T& y;
        T& z;

        // ---- Constructors ----
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

        // ---- Assigment ----
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



} // namespace april::utils