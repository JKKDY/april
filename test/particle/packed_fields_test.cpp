#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

#include "april/particle/source.hpp"
#include "april/particle/packed_access.hpp"
#include "april/base/types.hpp"

using namespace april;
using namespace april::env;

// Define a test mask (Standard Physics Fields)
static constexpr auto TestMask = Field::position | Field::velocity | Field::force | Field::mass;

class PackedParticleViewsTest : public testing::Test {
protected:
    static constexpr size_t Count = 16;
    
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<double> force_x, force_y, force_z;
    std::vector<double> mass;

    void SetUp() override {
        pos_x.resize(Count, 0.0); pos_y.resize(Count, 0.0); pos_z.resize(Count, 0.0);
        vel_x.resize(Count, 0.0); vel_y.resize(Count, 0.0); vel_z.resize(Count, 0.0);
        force_x.resize(Count, 0.0); force_y.resize(Count, 0.0); force_z.resize(Count, 0.0);
        mass.resize(Count, 0.0);

        // Particle i has position {i, i, i} and mass 1.0
        for(size_t i=0; i<Count; ++i) {
            const double val = static_cast<double>(i);
            pos_x[i] = val; pos_y[i] = val; pos_z[i] = val;
            vel_x[i] = 1.0; vel_y[i] = 0.0; vel_z[i] = 0.0; // Moving right
            mass[i] = 2.0;
        }
    }

    auto get_source() {
        ParticleSource<TestMask, NoUserData, false> src;
        
        // Setup SoA Pointers for Position
        src.position.x = pos_x.data();
        src.position.y = pos_y.data();
        src.position.z = pos_z.data();

        // Velocity
        src.velocity.x = vel_x.data();
        src.velocity.y = vel_y.data();
        src.velocity.z = vel_z.data();

        // Force
        src.force.x = force_x.data();
        src.force.y = force_y.data();
        src.force.z = force_z.data();

        // Mass
        src.mass = mass.data();

        return src;
    }
    
    // Helper: Verify data at specific index
    void ExpectParticle(const size_t index, const vec3& expected_pos, const vec3& expected_force) const {
        EXPECT_DOUBLE_EQ(pos_x[index], expected_pos.x);
        EXPECT_DOUBLE_EQ(pos_y[index], expected_pos.y);
        EXPECT_DOUBLE_EQ(pos_z[index], expected_pos.z);

        EXPECT_DOUBLE_EQ(force_x[index], expected_force.x);
        EXPECT_DOUBLE_EQ(force_y[index], expected_force.y);
        EXPECT_DOUBLE_EQ(force_z[index], expected_force.z);
    }
};


// Instantiation & Read Check
TEST_F(PackedParticleViewsTest, ReadValues) {
    const auto src = get_source();

    // Create the Packed Ref
    const PackedParticleRef<TestMask> p(src);

    // Read Position (Load)
    const pvec3 pos = p.position; // Implicit load of N particles

    // Verify Lanes
    // If SIMD width is 4, we expect {0,1,2,3} in the X lane
    const auto x_vals = pos.x.to_array();
    for(size_t i=0; i<x_vals.size(); ++i) {
        EXPECT_DOUBLE_EQ(x_vals[i], static_cast<double>(i));
    }
}

// Write Check (Broadcast)
TEST_F(PackedParticleViewsTest, WriteBroadcast) {
    const auto src = get_source();
    PackedParticleRef<TestMask> p(src);

    // Write constant force to all particles in the SIMD chunk
    p.force = pvec3(10.0, 20.0, 30.0);

    // Verify memory was updated for the first SIMD width
    constexpr size_t width = packedd::size();
    for(size_t i=0; i<width; ++i) {
        EXPECT_DOUBLE_EQ(force_x[i], 10.0);
        EXPECT_DOUBLE_EQ(force_y[i], 20.0);
        EXPECT_DOUBLE_EQ(force_z[i], 30.0);
    }

    // Verify index 'width' (next chunk) was NOT touched
    EXPECT_DOUBLE_EQ(force_x[width], 0.0);
}

// Physics Kernel Logic (Read-Modify-Write)
TEST_F(PackedParticleViewsTest, PhysicsUpdate) {
    const auto src = get_source();
    PackedParticleRef<TestMask> p(src);

    constexpr double dt = 0.1;

    // A simple Euler step: pos += vel * dt
    p.position += p.velocity * dt;

    // Verify Calculation
    constexpr size_t width = packedd::size();
    for(size_t i=0; i<width; ++i) {
        const double start_pos = static_cast<double>(i);
        const double expected = start_pos + 1.0 * dt; // vel.x is 1.0

        EXPECT_DOUBLE_EQ(pos_x[i], expected);

        // Y and Z shouldn't change (vel is 0)
        EXPECT_DOUBLE_EQ(pos_y[i], start_pos);
    }
}

// Force Kernel Logic (Interaction)
TEST_F(PackedParticleViewsTest, ForceKernel) {
    const auto src = get_source();
    PackedParticleRef<TestMask> p(src);

    // Simple drag force: F = -v * mass
    p.force = -p.velocity * p.mass;

    constexpr size_t width = packedd::size();
    for(size_t i=0; i<width; ++i) {
        // Expected: -1.0 * 2.0 = -2.0
        EXPECT_DOUBLE_EQ(force_x[i], -2.0);
        EXPECT_DOUBLE_EQ(force_y[i], 0.0);
        EXPECT_DOUBLE_EQ(force_z[i], 0.0);
    }
}



// 5. Const View Safety
TEST_F(PackedParticleViewsTest, ConstView) {
    const auto src = get_source();
    PackedParticleRef<TestMask> ref(src);

    // Convert to View
    const PackedParticleView<TestMask> view = ref.to_view();

    // Read is allowed
    const pvec3 v = view.velocity;
    EXPECT_DOUBLE_EQ(v.x.to_array()[0], 1.0);

    // Write should fail compile (uncomment to manually verify)
    // view.velocity = pvec3(0.0);
}