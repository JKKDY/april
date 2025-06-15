#pragma once
#include <cmath>
#include <string>
#include <format>
#include <concepts>

#include "april/utils/debug.h"


namespace april::utils {
    template <typename  T> requires std::integral<T> || std::floating_point<T>
    struct Vec3 {
        T x, y, z;

        Vec3() : x(0.0), y(0.0), z(0.0) {}
        explicit Vec3(const T v): x(v), y(v), z(v) {}
        Vec3(const T x, const T y, const T z) : x(x), y(y), z(z) {}

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
        Vec3& mul(const Vec3 & other) noexcept {
            x *= other.x;
            y *= other.y;
            z *= other.z;
            return *this;
        }

        // scalar product
        T operator*=(const Vec3 & other) const noexcept {
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
    };
} // namespace april::utils