#pragma once
#include <cmath>
#include <string>
#include <format>
#include <concepts>
#include <functional>

#include "april/utils/debug.h"


namespace april::utils {
    template <typename  T> requires std::integral<T> || std::floating_point<T>
    struct Vec3 {
        T x, y, z;

        Vec3() : x(static_cast<T>(0.0)), y(static_cast<T>(0.0)), z(static_cast<T>(0.0)) {}
        Vec3(const T x, const T y, const T z) : x(x), y(y), z(z) {}
        explicit Vec3(const T v): x(v), y(v), z(v) {}

        // unary minus
        Vec3 operator-() const noexcept {
            return Vec3{-x, -y, -z};
        }

        // Addition operator: returns a new vector that is the sum of this and another vector
        Vec3 operator+(const Vec3& other) const noexcept {
            return {x + other.x, y + other.y, z + other.z};
        }

        // Subtraction operator: returns a new vector that is the difference of this and another vector
        Vec3 operator-(const Vec3& other) const noexcept {
            return {x - other.x, y - other.y, z - other.z};
        }

        // Scalar multiplication operator: returns a new vector scaled by a constant
        Vec3 operator*(const T scalar) const noexcept {
            return {x * scalar, y * scalar, z * scalar};
        }

        // Scalar division: returns a new vector inversely scaled
        Vec3 operator/(const T scalar) const noexcept {
            return {x / scalar, y / scalar, z / scalar};
        }

        // Friend function for scalar multiplication with the scalar on the left
        friend Vec3 operator*(const T scalar, const Vec3& vec) noexcept {
            return vec * scalar;
        }

        // Compound assignment for addition
        Vec3& operator+=(const Vec3& other) noexcept {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        // Compound assignment for subtraction
        Vec3& operator-=(const Vec3& other) noexcept {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        // point-wise multiplication
        Vec3& operator*(const Vec3 & other) noexcept {
            x *= other.x;
            y *= other.y;
            z *= other.z;
            return *this;
        }

        Vec3 operator*(const Vec3 & other) const noexcept {
            return {x * other.x, y * other.y, z * other.z};
        }

        // point wise division
        Vec3& operator/(const Vec3 & other) noexcept {
            x /= other.x;
            y /= other.y;
            z /= other.z;
            return *this;
        }

        Vec3 operator/(const Vec3 & other) const noexcept {
            return {x / other.x, y / other.y, z / other.z};
        }

        Vec3 hadamard(const Vec3& other) { return *this * other; }
        Vec3 elementwise_div(const Vec3& other) { return *this / other; }


        // scalar product
        T dot(const Vec3 & other) const noexcept {
            return x * other.x + y * other.y + z * other.z;
        }

        // Compound assignment for scalar multiplication
        Vec3& operator*=(const T scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        // Access component by index: 0 for x, 1 for y, 2 for z.
        T operator[](const int index) const noexcept{
            AP_ASSERT(index >= 0 && index < 3, "Index out of bounds");
            switch (index) {
                case 0: return x;
                case 1: return y;
                case 2: return z;
            default: ;
            }
			return 0; // This line should never be reached due to the assertion above.
        }

        T & operator[](const int index) noexcept {
            AP_ASSERT(index >= 0 && index < 3, "Index out of bounds");
            switch (index) {
                case 0: return x;
                case 1: return y;
                case 2: return z;
            default: ;
            }
            return x; // This line should never be reached due to the assertion above.
        }

		[[nodiscard]] double norm_squared() const noexcept {
			return x * x + y * y + z * z;
		}

        [[nodiscard]] double inv_norm() const noexcept {
            return 1 / std::sqrt(norm_squared()); // compiler may optimize with fast inverse square root
        }

		[[nodiscard]] double norm() const noexcept {
			return sqrt(norm_squared());
		}

        bool operator==(const Vec3 & other) const noexcept {
            return x == other.x && y == other.y && z == other.z;
        }

        bool operator<=(const Vec3& other) const noexcept {
            return x <= other.x and y <= other.y and z <= other.z;
        }

        bool operator>=(const Vec3& other) const noexcept {
            return x >= other.x and y >= other.y and z >= other.z;
        }

        bool operator<(const Vec3& other) const noexcept {
            return x < other.x and y < other.y and z < other.z;
        }

        bool operator>(const Vec3& other) const noexcept {
            return x > other.x and y > other.y and z > other.z;
        }

        [[nodiscard]] std::string to_string() const {
            return std::format("{{{}, {}, {}}}", x, y, z);
        }

        T max() {
            return std::max(x, std::max(y, z));
        }

        T min() {
            return std::min(x, std::min(y, z));
        }

        static bool any(Vec3 v, std::function<bool(T)> c) {
            return c(v.x) || c(v.y) || c(v.z);
        }

        static bool all(Vec3 v, std::function<bool(T)> c) {
            return c(v.x) && c(v.y) && c(v.z);
        }
    };
} // namespace april::utils