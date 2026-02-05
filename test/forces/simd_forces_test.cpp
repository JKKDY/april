#include <gtest/gtest.h>
#include <vector>
#include <cmath>

#include "april/simd/packed.hpp"

#include "april/forces/lennard_jones.hpp"
#include "april/forces/harmonic.hpp"
#include "april/forces/gravity.hpp"
#include "april/forces/no_force.hpp"

using namespace april;

// 1. Define Test Types (SIMD Backends)
using PackedTypes = testing::Types<
    simd::Packed<double>
>;

template <typename T>
class ForceKernelTest : public testing::Test {
public:
    using Packed = T;
    using Scalar = T::value_type;
    using Vec3P  = math::Vec3<Packed>; // Packed Vector
    using Vec3S  = math::Vec3<Scalar>; // Scalar Vector

    // Helpers
    static constexpr size_t Width = Packed::size();

    // Helper: Compare SIMD result (Vec3P) vs Expected (vector<Vec3S>)
    void ExpectMatch(const Vec3P& simd_res, const std::vector<Vec3S>& expected, const char* name) {
        auto rx = simd_res.x.to_array();
        auto ry = simd_res.y.to_array();
        auto rz = simd_res.z.to_array();

        for (size_t i = 0; i < Width; ++i) {
            EXPECT_NEAR(rx[i], expected[i].x, 1e-13) << name << " X mismatch at lane " << i;
            EXPECT_NEAR(ry[i], expected[i].y, 1e-13) << name << " Y mismatch at lane " << i;
            EXPECT_NEAR(rz[i], expected[i].z, 1e-13) << name << " Z mismatch at lane " << i;
        }
    }
};

TYPED_TEST_SUITE(ForceKernelTest, PackedTypes);

// ============================================================================
// 1. LENNARD-JONES (12-6 Potential)
// ============================================================================
TYPED_TEST(ForceKernelTest, LennardJones) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Vec3P  = TestFixture::Vec3P;
    using Vec3S  = TestFixture::Vec3S;

    force::LennardJones lj(1.0, 1.0); // Epsilon=1, Sigma=1

    // Setup Inputs
    std::vector<Vec3S> inputs(TestFixture::Width);
    std::vector<Vec3S> expected(TestFixture::Width);

    std::vector<Scalar> ix(TestFixture::Width), iy(TestFixture::Width), iz(TestFixture::Width);

    for(size_t i=0; i<TestFixture::Width; ++i) {
        // Vary distance from 0.9 to 1.5
        double dist = 0.9 + (0.6 * i) / TestFixture::Width;
        inputs[i] = Vec3S(dist, 0.0, 0.0);

        // Scalar Calc
        expected[i] = lj.eval(0, 0, inputs[i]);

        // Prep SIMD Data
        ix[i] = inputs[i].x; iy[i] = inputs[i].y; iz[i] = inputs[i].z;
    }

    // Load SIMD
    Vec3P r_packed(Packed::load(ix.data()), Packed::load(iy.data()), Packed::load(iz.data()));
    // Vec3P r = {-r_packed.x, -r_packed.y, -r_packed.z};
    Vec3P res = lj.eval(0, 0, r_packed);
    // auto x = -l;
    // Run SIMD Kernel

    this->ExpectMatch(res, expected, "LennardJones");
}

// ============================================================================
// 2. HARMONIC (Spring Force)
// ============================================================================
TYPED_TEST(ForceKernelTest, HarmonicSpring) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Vec3P  = TestFixture::Vec3P;
    using Vec3S  = TestFixture::Vec3S;

    constexpr double k = 100.0;
    constexpr double r0 = 2.0;
    force::Harmonic spring(k, r0);

    std::vector<Vec3S> inputs(TestFixture::Width);
    std::vector<Vec3S> expected(TestFixture::Width);
    std::vector<Scalar> ix(TestFixture::Width), iy(TestFixture::Width), iz(TestFixture::Width);

    for(size_t i=0; i<TestFixture::Width; ++i) {
        // Distances: Some compressed (<2.0), some stretched (>2.0)
        const double dist = 1.5 + (1.0 * i) / TestFixture::Width; // 1.5 to 2.5

        // Use diagonals to test norm() logic
        double axis = dist / std::sqrt(2.0); // x=axis, y=axis -> norm=dist
        inputs[i] = Vec3S(axis, axis, 0.0);

        // Scalar Calc
        expected[i] = spring.eval(0, 0, inputs[i]);

        ix[i] = inputs[i].x; iy[i] = inputs[i].y; iz[i] = inputs[i].z;
    }

    Vec3P r_packed(Packed::load(ix.data()), Packed::load(iy.data()), Packed::load(iz.data()));

    // Run SIMD Kernel
    Vec3P res = spring.eval(0, 0, r_packed);

    this->ExpectMatch(res, expected, "Harmonic");
}

// ============================================================================
// 3. GRAVITY (Newtonian 1/r^2)
// ============================================================================
TYPED_TEST(ForceKernelTest, NewtonianGravity) {
    using Packed = TestFixture::Packed;
    using Scalar = TestFixture::Scalar;
    using Vec3P  = TestFixture::Vec3P;
    using Vec3S  = TestFixture::Vec3S;

    // Use a distinct constant to verify math (G=9.81)
    force::Gravity grav(9.81);

    // 1. Define Mock Particles
    // We need different types for the Scalar loop vs the SIMD call.

    struct MockScalarParticle { Scalar mass; };
    struct MockPackedParticle { Packed mass; };

    // Create particles with specific masses to verify the mass multiplication logic
    MockScalarParticle p1_s{100.0};
    MockScalarParticle p2_s{5.0};

    // For SIMD, we "broadcast" the scalar mass into the Packed type
    MockPackedParticle p1_p{Packed(100.0)};
    MockPackedParticle p2_p{Packed(5.0)};

    // 2. Prepare Data
    std::vector<Vec3S> inputs(TestFixture::Width);
    std::vector<Vec3S> expected(TestFixture::Width);
    std::vector<Scalar> ix(TestFixture::Width), iy(TestFixture::Width), iz(TestFixture::Width);

    for(size_t i=0; i<TestFixture::Width; ++i) {
        double dist = 2.0 + static_cast<double>(i);
        inputs[i] = Vec3S(0.0, dist, 0.0);

        // Scalar Calc: Pass Scalar Mock Particles
        expected[i] = grav.eval(p1_s, p2_s, inputs[i]);

        ix[i] = inputs[i].x; iy[i] = inputs[i].y; iz[i] = inputs[i].z;
    }

    // 3. Run SIMD Kernel
    Vec3P r_packed(Packed::load(ix.data()), Packed::load(iy.data()), Packed::load(iz.data()));

    // SIMD Calc: Pass Packed Mock Particles
    // This ensures (p1.mass * p2.mass) is calculated as a SIMD vector multiplication
    Vec3P res = grav.eval(p1_p, p2_p, r_packed);

    this->ExpectMatch(res, expected, "NewtonianGravity");
}


TYPED_TEST(ForceKernelTest, NoForce) {
    using Packed = TestFixture::Packed;
    using Vec3P  = TestFixture::Vec3P;

    force::NoForce no_force;

    // Inputs: (Values don't matter, but shouldn't crash)
    // We use huge numbers just to be sure it doesn't try to normalize them
    Vec3P r_packed(Packed(1e20), Packed(1e20), Packed(1e20));

    Vec3P res = no_force.eval(0, 0, r_packed);

    auto rx = res.x.to_array();
    auto ry = res.y.to_array();
    auto rz = res.z.to_array();

    for(size_t i=0; i<TestFixture::Width; ++i) {
        EXPECT_DOUBLE_EQ(rx[i], 0.0);
        EXPECT_DOUBLE_EQ(ry[i], 0.0);
        EXPECT_DOUBLE_EQ(rz[i], 0.0);
    }
}