#pragma once
#include <math.h>

#include "Debug.h"


namespace april::utils {
    struct Vec3 {
        double x, y, z;

        Vec3() : x(0.0), y(0.0), z(0.0) {}
        Vec3(double v): x(v), y(v), z(v) {}
        Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

        // Addition operator: returns a new vector that is the sum of this and another vector
        Vec3 operator+(const Vec3& other) const {
            return Vec3(x + other.x, y + other.y, z + other.z);
        }

        // Subtraction operator: returns a new vector that is the difference of this and another vector
        Vec3 operator-(const Vec3& other) const {
            return Vec3(x - other.x, y - other.y, z - other.z);
        }

        // Scalar multiplication operator: returns a new vector scaled by a constant
        inline Vec3 operator*(double scalar) const {
            return Vec3(x * scalar, y * scalar, z * scalar);
        }

        // Friend function for scalar multiplication with the scalar on the left
        friend inline Vec3 operator*(double scalar, const Vec3& vec) {
            return vec * scalar;
        }

        // Compound assignment for addition
        Vec3& operator+=(const Vec3& other) {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        // Compound assignment for subtraction
        Vec3& operator-=(const Vec3& other) {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        // pointwise multiplication multiplication
        Vec3& mul(const Vec3 & other) {
            x *= other.x;
            y *= other.y;
            z *= other.z;
            return *this;
        }

        // scalar product
        double operator*=(const Vec3 & other) {
            return x * other.x + y * other.y + z * other.z;
        }

        // Compound assignment for scalar multiplication
        Vec3& operator*=(double scalar) {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        // Access component by index
        // index: 0 for x, 1 for y, 2 for z.
        double operator[](int index) const {
            AP_ASSERT(index >= 0 && index < 3, "Index out of bounds");
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            }
        }

        double & operator[](int index) {
            AP_ASSERT(index >= 0 && index < 3, "Index out of bounds");
            switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            }
        }

		double norm_squared() const {
			return x * x + y * y + z * z;
		}

		double norm() const {
			return sqrt(norm_squared());
		}

        bool operator==(const Vec3 & other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };
} // namespace april::utils