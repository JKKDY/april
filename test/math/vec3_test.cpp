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

    Vec3T v(1.0, 2.0, 3.0);

    this->ExpectEq(v.x, 1.0);
    this->ExpectEq(v.y, 2.0);
    this->ExpectEq(v.z, 3.0);
}

// Basic Arithmetic (Vector-Vector)
TYPED_TEST(Vec3Test, Arithmetic) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;

    Vec3T a(1.0, 2.0, 3.0);
    Vec3T b(4.0, 5.0, 6.0);

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

    Vec3T v(1.0, 2.0, 3.0);

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
    Vec3T v(0.0, 3.0, 4.0);

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

    Vec3T v(10.0, 10.0, 10.0);
    Vec3T u(1.0, 2.0, 3.0);

    v += u;
    this->ExpectEq(v.x, 11.0);
    this->ExpectEq(v.y, 12.0);

    v *= 2.0;
    this->ExpectEq(v.x, 22.0);
}

// Mixed Vector Arithmetic (Vec3<Wide> vs Vec3<double>)
TYPED_TEST(Vec3Test, MixedVectorArithmetic) {
    using T = TypeParam;                 // The Testing Type (likely Wide<double>)
    using Vec3T = april::math::Vec3<T>;  // The Target Vector (Vec3<Wide>)
    using Vec3S = april::math::Vec3<double>; // The Source Scalar Vector

    // Setup
    Vec3S scalar_vec(1.0, 2.0, 3.0);       // "Global" vector (e.g. Gravity, Wind)
    Vec3T wide_vec(10.0, 20.0, 30.0); // "Particle" vector

    // Construction from Scalar Vector (Broadcast)
    // "Vec3<Wide> w = v_scalar;"
    Vec3T converted(scalar_vec);
    this->ExpectEq(converted.x, 1.0);
    this->ExpectEq(converted.y, 2.0);
    this->ExpectEq(converted.z, 3.0);

    // Addition (Wide + ScalarVec)
    Vec3T sum = wide_vec + scalar_vec;
    this->ExpectEq(sum.x, 11.0); // 10 + 1
    this->ExpectEq(sum.y, 22.0); // 20 + 2
    this->ExpectEq(sum.z, 33.0); // 30 + 3

    // Subtraction (Wide - ScalarVec)
    Vec3T diff = wide_vec - scalar_vec;
    this->ExpectEq(diff.x, 9.0);
    this->ExpectEq(diff.y, 18.0);

    // Hadamard Product (Wide * ScalarVec)
    Vec3T prod = wide_vec * scalar_vec;
    this->ExpectEq(prod.x, 10.0); // 10 * 1
    this->ExpectEq(prod.y, 40.0); // 20 * 2
    this->ExpectEq(prod.z, 90.0); // 30 * 3

    // Compound Assignment (Wide += ScalarVec)
    wide_vec += scalar_vec;
    this->ExpectEq(wide_vec.x, 11.0);
    this->ExpectEq(wide_vec.y, 22.0);
}

// Complex Expression "Mish Mash"
// Verifies that temporary objects and mixed types survive a long equation
TYPED_TEST(Vec3Test, ArithmeticMishMash) {
    using T = TypeParam;
    using Vec3T = math::Vec3<T>;
    using Vec3S = math::Vec3<double>;

    Vec3T pos(10.0, 10.0, 10.0);
    Vec3T old_pos(9.0, 10.0, 11.0);
    Vec3S gravity(0.0, -10.0, 0.0);
    constexpr double dt = 0.1;
    constexpr double damping = 0.99;

    // A standard Verlet-integration-like expression:
    // next_pos = pos + (pos - old_pos) * damping + gravity * (dt * dt)

    // Breakdown:
    // 1. (pos - old_pos): Vec3<Wide> result -> {1, 0, -1}
    // 2. * damping:       Vec3<Wide> result -> {0.99, 0, -0.99}
    // 3. gravity * dt*dt: Vec3<Scalar> * double -> Vec3<Scalar> {0, -0.1, 0}
    // 4. Vec3<Wide> + Vec3<Scalar>: Broadcast add

    Vec3T next_pos = pos + (pos - old_pos) * damping + gravity * (dt * dt);

    this->ExpectEq(next_pos.x, 10.99);        // 10 + 0.99 + 0
    this->ExpectEq(next_pos.y, 9.9);          // 10 + 0    - 0.1
    this->ExpectEq(next_pos.z, 9.01);         // 10 - 0.99 + 0
}