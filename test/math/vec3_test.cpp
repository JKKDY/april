#include <gtest/gtest.h>
#include <cmath>

#include "april/base/types.hpp"
#include "april/simd/backend_xsimd.hpp" 
#include "april/simd/backend_std_simd.hpp"

using namespace april;

// We test:
// - double
// - Wide<double> (SIMD physics)
using Vec3Types = testing::Types<
    double, 
    simd::internal::xsimd::Wide<double>,
    simd::internal::std_simd::Wide<double>
>;

template <typename T>
class Vec3Test : public testing::Test {
public:
    // Helper to abstract the "Verification" step.
    // If T is double, we check equality directly.
    // If T is Wide, we check all lanes.
    void ExpectEq(const T& actual, double expected_val, double tolerance = 1e-15) {
        if constexpr (std::is_same_v<T, double>) {
            EXPECT_DOUBLE_EQ(actual, expected_val);
        } else {
            // Use our new to_array helper
            auto arr = actual.to_array();
            for (size_t i = 0; i < arr.size(); ++i) {
                EXPECT_NEAR(arr[i], expected_val, tolerance)
                    << "Mismatch at SIMD lane " << i;
            }
        }
    }

    // Helper to construct a value
    // If T is Wide, it broadcasts 'val'.
    T Val(double val) { return T(val); }
};

TYPED_TEST_SUITE(Vec3Test, Vec3Types);

// Construction & Component Access
TYPED_TEST(Vec3Test, Construction) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    Vec3T v(this->Val(1.0), this->Val(2.0), this->Val(3.0));

    this->ExpectEq(v.x, 1.0);
    this->ExpectEq(v.y, 2.0);
    this->ExpectEq(v.z, 3.0);
}

// Basic Arithmetic (Vector-Vector)
TYPED_TEST(Vec3Test, Arithmetic) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    Vec3T a(this->Val(1.0), this->Val(2.0), this->Val(3.0));
    Vec3T b(this->Val(4.0), this->Val(5.0), this->Val(6.0));

    // Addition
    Vec3T sum = a + b;
    this->ExpectEq(sum.x, 5.0);
    this->ExpectEq(sum.y, 7.0);
    this->ExpectEq(sum.z, 9.0);

    // Subtraction
    Vec3T diff = b - a;
    this->ExpectEq(diff.x, 3.0);
    
    // Element-wise Multiplication (Hadamard)
    Vec3T prod = a * b; // 1*4, 2*5, 3*6
    this->ExpectEq(prod.x, 4.0);
    this->ExpectEq(prod.y, 10.0);
    this->ExpectEq(prod.z, 18.0);
}

// Scalar Operations (Mixed Arithmetic)
TYPED_TEST(Vec3Test, ScalarOps) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    Vec3T v(this->Val(1.0), this->Val(2.0), this->Val(3.0));

    // Vec * Scalar
    Vec3T scaled = v * 2.0; 
    this->ExpectEq(scaled.x, 2.0);
    this->ExpectEq(scaled.y, 4.0);
    this->ExpectEq(scaled.z, 6.0);

    // Scalar * Vec (Friend operator)
    Vec3T scaled_left = 3.0 * v;
    this->ExpectEq(scaled_left.x, 3.0);
    
    // Division
    Vec3T div = v / 2.0;
    this->ExpectEq(div.x, 0.5);
}

// Geometry Ops (Dot, Norm)
TYPED_TEST(Vec3Test, Geometry) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    // 3-4-5 Triangle vector
    Vec3T v(this->Val(0.0), this->Val(3.0), this->Val(4.0));

    // Norm Squared
    T n2 = v.norm_squared();
    this->ExpectEq(n2, 25.0); // 0 + 9 + 16

    // Norm (Requires sqrt ADL to work)
    T n = v.norm();
    this->ExpectEq(n, 5.0);

    // Inverse Norm (Requires rsqrt ADL/Fallback)
    // We allow a looser tolerance for the hardware approximation
    T inv = v.inv_norm();
    this->ExpectEq(inv, 0.2, 1e-3);
}

// Compound Assignment
TYPED_TEST(Vec3Test, CompoundAssignment) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    Vec3T v(this->Val(10.0), this->Val(10.0), this->Val(10.0));
    Vec3T u(this->Val(1.0), this->Val(2.0), this->Val(3.0));

    v += u;
    this->ExpectEq(v.x, 11.0);
    this->ExpectEq(v.y, 12.0);

    v *= 2.0;
    this->ExpectEq(v.x, 22.0);
}