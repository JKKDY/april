#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <cmath>

#include "april/base/types.hpp"
#include "april/simd/backend_xsimd.hpp"
#include "april/simd/backend_std_simd.hpp"

using namespace april;

// The SIMD Types to test
using SimdTypes = testing::Types<
    simd::internal::xsimd::Packed<double>,
    simd::internal::std_simd::Packed<double>
>;

template <typename T>
class SimdProxyTest : public testing::Test {
public:
    using Scalar = T::value_type;
    using Vec3T  = math::Vec3<T>;           // Vec3<Packed>
    using Proxy  = math::Vec3Proxy<T>;      // The SoA Proxy
    using Vec3S  = math::Vec3<Scalar>;      // Scalar Vec3 (gravity, etc.)

    // Backing Memory (SoA Layout)
    // We resize these to be large enough for the widest SIMD vector (e.g., AVX512 = 8 doubles)
    std::vector<Scalar> x_buf, y_buf, z_buf;

    void SetUp() override {
        // Ensure strictly aligned memory or just enough space.
        // Vectors are usually aligned enough for load_unaligned.
        size_t size = std::max<size_t>(16, T::size());
        x_buf.resize(size, 0.0);
        y_buf.resize(size, 0.0);
        z_buf.resize(size, 0.0);
    }

    // Create a Proxy pointing to the start of the buffers
    Proxy MakeProxy() {
        return Proxy(x_buf.data(), y_buf.data(), z_buf.data());
    }

    // Helper: Verify that EVERY lane in the SIMD width matches the expected value
    void ExpectAllLanes(Scalar x, Scalar y, Scalar z) {
        for(size_t i=0; i < T::size(); ++i) {
            EXPECT_DOUBLE_EQ(x_buf[i], x) << "Mismatch in X at lane " << i;
            EXPECT_DOUBLE_EQ(y_buf[i], y) << "Mismatch in Y at lane " << i;
            EXPECT_DOUBLE_EQ(z_buf[i], z) << "Mismatch in Z at lane " << i;
        }
    }
};

TYPED_TEST_SUITE(SimdProxyTest, SimdTypes);

// 1. Write Value to Proxy -> Memory
// Verifies broadcasting: p.pos = {1, 2, 3} writes 1,2,3 to ALL lanes
TYPED_TEST(SimdProxyTest, BroadcastWrite) {
    auto p = this->MakeProxy();
    using Vec3T = TestFixture::Vec3T;

    // Implicit broadcast of scalars to Packed
    p = Vec3T(1.0, 2.0, 3.0);

    this->ExpectAllLanes(1.0, 2.0, 3.0);
}

// 2. Read Proxy -> Value
// Verifies Loading from memory
TYPED_TEST(SimdProxyTest, ReadFromMemory) {
    auto p = this->MakeProxy();

    // Setup Memory manually
    std::fill(this->x_buf.begin(), this->x_buf.end(), 5.0);
    std::fill(this->y_buf.begin(), this->y_buf.end(), 6.0);
    std::fill(this->z_buf.begin(), this->z_buf.end(), 7.0);

    // Read
    typename TestFixture::Vec3T v = p;

    // Verify the read Packed value (by storing it back to a temp buffer or checking lanes)
    // Here we just use the proxy to write it back to check logic
    // But better to check the value v directly.
    // Assuming Packed has to_array or similar, or we just trust the math tests.
    // Let's reuse the memory check by writing it back to zeroed memory.

    std::fill(this->x_buf.begin(), this->x_buf.end(), 0.0);
    std::fill(this->y_buf.begin(), this->y_buf.end(), 0.0);
    std::fill(this->z_buf.begin(), this->z_buf.end(), 0.0);

    p = v; // Write back
    this->ExpectAllLanes(5.0, 6.0, 7.0);
}

// 3. Component-wise Access
// p.pos.y = 99.0
TYPED_TEST(SimdProxyTest, ComponentAccess) {
    auto p = this->MakeProxy();

    // Assign to just the Y component (SimdRef assignment)
    p.y = 99.0;

    // Verify X and Z are 0, Y is 99
    this->ExpectAllLanes(0.0, 99.0, 0.0);
}

// 4. Arithmetic with Scalar Vector
// p.pos += Vec3d(1, 2, 3)
TYPED_TEST(SimdProxyTest, AddScalarVector) {
    auto p = this->MakeProxy();
    using Vec3S = TestFixture::Vec3S;

    // Init memory to 10
    p = typename TestFixture::Vec3T(10.0, 10.0, 10.0);

    Vec3S offset(1.0, 2.0, 3.0);

    p += offset; // Should broadcast offset adds to all lanes

    this->ExpectAllLanes(11.0, 12.0, 13.0);
}

// 5. Arithmetic with Pure Scalar
// p.pos *= 2.0
TYPED_TEST(SimdProxyTest, ScaleByScalar) {
    auto p = this->MakeProxy();

    // Init memory to {1, 2, 3}
    p = typename TestFixture::Vec3T(1.0, 2.0, 3.0);

    p *= 2.0;

    this->ExpectAllLanes(2.0, 4.0, 6.0);
}

// 6. Physics Integration (The Big One)
// p.pos += v*dt + g*0.5*dt*dt
TYPED_TEST(SimdProxyTest, PhysicsKernel) {
    auto p = this->MakeProxy();
    using Vec3T = TestFixture::Vec3T;
    using Vec3S = TestFixture::Vec3S;

    // Init Position: {10, 10, 10}
    p = Vec3T(10.0, 10.0, 10.0);

    // Velocity: {1, 0, 0} (Particle moving right)
    Vec3T velocity(1.0, 0.0, 0.0);

    // Gravity: {0, -10, 0} (Global constant)
    Vec3S gravity(0.0, -10.0, 0.0);

    double dt = 0.1;

    // Update
    p += velocity * dt + gravity * (0.5 * dt * dt);

    // Expected Results:
    // X: 10 + 1*0.1 + 0 = 10.1
    // Y: 10 + 0 - 10*0.5*0.01 = 9.95
    // Z: 10
    this->ExpectAllLanes(10.1, 9.95, 10.0);
}