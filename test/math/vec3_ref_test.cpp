#include <gtest/gtest.h>
#include "april/base/types.hpp"

using namespace april;

class ScalarProxyTest : public testing::Test {
public:
    // 1. Backing Memory (Simulating a particle's data)
    double x_mem = 0.0;
    double y_mem = 0.0;
    double z_mem = 0.0;

    // 2. The Proxy Type
    using Proxy = math::Vec3Proxy<double>;
    using Vec3  = math::Vec3<double>;

    // 3. Helper to create a proxy binding to our memory
    Proxy GetProxy() {
        return Proxy(x_mem, y_mem, z_mem);
    }

    void SetUp() override {
        // Reset memory before each test
        x_mem = 0.0; y_mem = 0.0; z_mem = 0.0;
    }
};

// 1. Basic Read/Write Verification
TEST_F(ScalarProxyTest, ReadWriteMemory) {
    auto p = GetProxy();

    // Write via Proxy -> Check Memory
    p = Vec3(10.0, 20.0, 30.0);
    EXPECT_DOUBLE_EQ(x_mem, 10.0);
    EXPECT_DOUBLE_EQ(y_mem, 20.0);
    EXPECT_DOUBLE_EQ(z_mem, 30.0);

    // Write Memory -> Check Proxy Read
    y_mem = 55.0;
    EXPECT_DOUBLE_EQ(p.y, 55.0);
}

// 2. Component Access
TEST_F(ScalarProxyTest, ComponentAccess) {
    auto p = GetProxy();

    // "p.pos.x = 5"
    p.x = 5.0;
    EXPECT_DOUBLE_EQ(x_mem, 5.0);

    // Ensure others are untouched
    EXPECT_DOUBLE_EQ(y_mem, 0.0);
}

// 3. Basic Arithmetic (Proxy += Vec3)
TEST_F(ScalarProxyTest, VectorArithmetic) {
    auto p = GetProxy();
    p = Vec3(10.0, 10.0, 10.0);

    Vec3 velocity(1.0, 2.0, 3.0);

    // p.pos += velocity
    p += velocity;

    EXPECT_DOUBLE_EQ(x_mem, 11.0);
    EXPECT_DOUBLE_EQ(y_mem, 12.0);
    EXPECT_DOUBLE_EQ(z_mem, 13.0);
}

// 4. Scalar Arithmetic (Proxy *= double)
TEST_F(ScalarProxyTest, ScalarArithmetic) {
    auto p = GetProxy();
    p = Vec3(2.0, 4.0, 8.0);

    // p.pos *= 0.5
    p *= 0.5;

    EXPECT_DOUBLE_EQ(x_mem, 1.0);
    EXPECT_DOUBLE_EQ(y_mem, 2.0);
    EXPECT_DOUBLE_EQ(z_mem, 4.0);
}

// 5. Proxy to Proxy Assignment
TEST_F(ScalarProxyTest, ProxyToProxy) {
    // Source: Particle A {1, 2, 3}
    double ax = 1.0, ay = 2.0, az = 3.0;
    Proxy pA(ax, ay, az);

    // Dest: Particle B (Our member vars)
    auto pB = GetProxy();

    // "pB.pos = pA.pos"
    pB = pA;

    EXPECT_DOUBLE_EQ(x_mem, 1.0);
    EXPECT_DOUBLE_EQ(y_mem, 2.0);
    EXPECT_DOUBLE_EQ(z_mem, 3.0);
}

// 6. Physics Expression Integration
TEST_F(ScalarProxyTest, PhysicsExpression) {
    auto p = GetProxy();

    // Initial State
    p = Vec3(0.0, 10.0, 0.0);

    Vec3 velocity(1.0, 0.0, 0.0);
    Vec3 gravity(0.0, -9.81, 0.0);
    double dt = 0.1;

    // "p.pos += v*dt + g*0.5*dt^2"
    p += velocity * dt + gravity * (0.5 * dt * dt);

    // Expected X: 0 + 1*0.1 = 0.1
    EXPECT_DOUBLE_EQ(x_mem, 0.1);

    // Expected Y: 10 + 0 - 9.81*0.5*0.01 = 10 - 0.04905 = 9.95095
    EXPECT_DOUBLE_EQ(y_mem, 9.95095);
}