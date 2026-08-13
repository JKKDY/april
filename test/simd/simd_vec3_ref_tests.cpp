#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

#include "april/base/types.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/masked_packed.hpp"

using namespace april;




using SimdTypes = testing::Types<
    april::simd::Packed<float, april::simd::float_width>,
    april::simd::Packed<double, april::simd::double_width>
>;


template<typename T>
class SimdProxyTest : public testing::Test {
public:
    using Packed = T;
    using Scalar = Packed::value_type;

    using Location = simd::ContiguousLocation<Scalar>;
    using Ref = simd::PackedRef<Location>;

    using Vec3T = math::Vec3<Packed>;
    using Proxy = math::Vec3Proxy<Ref>;
    using Vec3S = math::Vec3<Scalar>;

    std::vector<Scalar> x_buf, y_buf, z_buf;

    void SetUp() override {
        const size_t size = std::max<size_t>(16, Packed::size());

        x_buf.resize(size, Scalar{0});
        y_buf.resize(size, Scalar{0});
        z_buf.resize(size, Scalar{0});
    }

    Proxy MakeProxy() {
        return Proxy{
            Ref{Location{x_buf.data()}},
            Ref{Location{y_buf.data()}},
            Ref{Location{z_buf.data()}}
        };
    }

    void ExpectAllLanes(Scalar x, Scalar y, Scalar z) {
        for (size_t i = 0; i < Packed::size(); ++i) {
            EXPECT_DOUBLE_EQ(x_buf[i], x)
                << "Mismatch in X at lane " << i;

            EXPECT_DOUBLE_EQ(y_buf[i], y)
                << "Mismatch in Y at lane " << i;

            EXPECT_DOUBLE_EQ(z_buf[i], z)
                << "Mismatch in Z at lane " << i;
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
    std::fill(this->x_buf.begin(), this->x_buf.end(), static_cast<TestFixture::Scalar>(5.0));
    std::fill(this->y_buf.begin(), this->y_buf.end(), static_cast<TestFixture::Scalar>(6.0));
    std::fill(this->z_buf.begin(), this->z_buf.end(), static_cast<TestFixture::Scalar>(7.0));

    // Read
    typename TestFixture::Vec3T v = p;

    // Verify the read Packed value (by storing it back to a temp buffer or checking lanes)
    // Here we just use the proxy to write it back to check logic
    // But better to check the value v directly.
    // Assuming Packed has to_array or similar, or we just trust the math tests.
    // Let's reuse the memory check by writing it back to zeroed memory.

    std::fill(this->x_buf.begin(), this->x_buf.end(), static_cast<TestFixture::Scalar>(5.0));
    std::fill(this->y_buf.begin(), this->y_buf.end(), static_cast<TestFixture::Scalar>(6.0));
    std::fill(this->z_buf.begin(), this->z_buf.end(), static_cast<TestFixture::Scalar>(7.0));

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
    using Scalar = TestFixture::Scalar;

    // Init Position: {10, 10, 10}
    p = Vec3T(10.0, 10.0, 10.0);

    // Velocity: {1, 0, 0} (Particle moving right)
    Vec3T velocity(1.0, 0.0, 0.0);

    // Gravity: {0, -10, 0} (Global constant)
    Vec3S gravity(
        static_cast<Scalar>(0.0),
        static_cast<Scalar>(-9.81),
        static_cast<Scalar>(0.0)
    );

    double dt = 0.1;

    // Update
    p += velocity * dt + gravity * (0.5 * dt * dt);

    // Expected Results:
    // X: 10 + 1*0.1 + 0 = 10.1
    // Y: 10 + 0 - 10*0.5*0.01 = 9.95
    // Z: 10
    this->ExpectAllLanes(
        static_cast<Scalar>(10.1),
        static_cast<Scalar>(9.95095),
        static_cast<Scalar>(10.0)
    );
}



// Vec3Proxy<Packed> += Vec3<MaskedPacked<Packed>>
TYPED_TEST(SimdProxyTest, AddMaskedPackedVectorToProxy) {
    using Packed = TypeParam;
    using Mask = typename Packed::mask_type;
    using Masked = simd::MaskedPacked<Packed, Mask>;
    using MaskedVec3 = math::Vec3<Masked>;

    auto p = this->MakeProxy();
    p = typename TestFixture::Vec3T(10.0, 20.0, 30.0);

    const Mask mask(true);

    MaskedVec3 delta(
        Masked(Packed(1.0), mask),
        Masked(Packed(2.0), mask),
        Masked(Packed(3.0), mask)
    );

    p += delta;

    this->ExpectAllLanes(11.0, 22.0, 33.0);
}


TYPED_TEST(SimdProxyTest, AddLaneVaryingMaskedPackedVectorToProxy) {
    using Packed = TypeParam;
    using Scalar = typename TestFixture::Scalar;
    using Mask = typename Packed::mask_type;
    using Masked = simd::MaskedPacked<Packed, Mask>;
    using MaskedVec3 = math::Vec3<Masked>;

    auto p = this->MakeProxy();
    p = typename TestFixture::Vec3T(10.0, 20.0, 30.0);

    std::vector<Scalar> x_values(Packed::size());
    std::vector<Scalar> y_values(Packed::size());
    std::vector<Scalar> z_values(Packed::size());

    for (size_t i = 0; i < Packed::size(); ++i) {
        x_values[i] = static_cast<Scalar>(i + 1);
        y_values[i] = static_cast<Scalar>(i + 2);
        z_values[i] = static_cast<Scalar>(i + 3);
    }

    const Mask mask(true);

    MaskedVec3 delta(
        Masked(Packed::load_unaligned(x_values.data()), mask),
        Masked(Packed::load_unaligned(y_values.data()), mask),
        Masked(Packed::load_unaligned(z_values.data()), mask)
    );

    p += delta;

    for (size_t i = 0; i < Packed::size(); ++i) {
        EXPECT_DOUBLE_EQ(this->x_buf[i], 10.0 + x_values[i])
            << "Mismatch in X at lane " << i;

        EXPECT_DOUBLE_EQ(this->y_buf[i], 20.0 + y_values[i])
            << "Mismatch in Y at lane " << i;

        EXPECT_DOUBLE_EQ(this->z_buf[i], 30.0 + z_values[i])
            << "Mismatch in Z at lane " << i;
    }
}










