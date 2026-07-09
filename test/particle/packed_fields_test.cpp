#include <gtest/gtest.h>
#include <vector>
#include <algorithm>

#include "../../include/april/particle/access/source.hpp"
#include "../../include/april/particle/access/packed_access.hpp"
#include "april/base/types.hpp"

using namespace april;
using namespace april::particle;
using namespace april::particle::internal;

// Define a test mask (Standard Physics Fields)
static constexpr auto TestMask = ParticleField::position | ParticleField::velocity | ParticleField::force | ParticleField::mass;

using PackedParticle = PackedParticleRef<TestMask, TestMask, NoParticleAttributes>;

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
        internal::ParticleSource<TestMask, TestMask, NoParticleAttributes> src;

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
    const PackedParticle p(src);

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
    PackedParticle p(src);

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
   PackedParticle p(src);

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
    PackedParticle p(src);

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
template <typename T>
concept VelocityAssignable =
    requires(T t) {
    t.velocity = pvec3(0.0);
    };
TEST_F(PackedParticleViewsTest, ConstView) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask,NoParticleAttributes> ref(src);

    // Convert to View
    const PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes> view = ref.to_view();

    // Read is allowed
    const pvec3 v = view.velocity;
    EXPECT_DOUBLE_EQ(v.x.to_array()[0], 1.0);

    // Write should fail compile on view, should succeed on ref
    // TODO
    // static_assert(!VelocityAssignable<PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes>>);
    static_assert( VelocityAssignable<PackedParticleRef<TestMask, TestMask, NoParticleAttributes>>);
}



// 6. Buffer Load & Isolation
// Verifies that 'load_buffer' reads correctly across ALL lanes,
// and that modifying the registers DOES NOT touch memory until commit.
TEST_F(PackedParticleViewsTest, BufferIsolation) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask, NoParticleAttributes> ref(src);

    // 1. Load into Registers
    auto buffer = ref.load_buffer();

    // Verify initial load (Memory was set to Vel={1.0, 0.0, 0.0})
    const auto buf_vel_x = buffer.velocity.x.to_array();
    constexpr size_t Width = packedd::size();

    for(size_t i = 0; i < Width; ++i) {
        EXPECT_DOUBLE_EQ(buf_vel_x[i], 1.0) << "Initial load mismatch at lane " << i;
    }

    // 2. Modify Register Only (Set all lanes to 99.0)
    buffer.velocity = pvec3(99.0, 99.0, 99.0);

    // 3. Verify Memory is UNTOUCHED
    // We explicitly read from the reference (proxy) to ensure it hits memory
    const auto ref_vel_x = ref.velocity.x.to_array();

    for(size_t i = 0; i < Width; ++i) {
        // Memory should still be 1.0
        EXPECT_DOUBLE_EQ(ref_vel_x[i], 1.0) << "Memory corruption at lane " << i;
        // Backing store check
        EXPECT_DOUBLE_EQ(vel_x[i], 1.0);
    }

    // 4. Verify Register IS Modified
    const auto new_buf_vel_x = buffer.velocity.x.to_array();
    for(size_t i = 0; i < Width; ++i) {
        EXPECT_DOUBLE_EQ(new_buf_vel_x[i], 99.0) << "Register update failed at lane " << i;
    }
}

// 7. Buffer Rotation
// Verifies "rotate_left" shifts elements correctly across the whole vector.
TEST_F(PackedParticleViewsTest, BufferRotation) {
    const auto src = get_source();
    const PackedParticleRef<TestMask, TestMask, NoParticleAttributes> ref(src);

    // Memory setup: Position X is {0, 1, 2, 3 ...}
    auto buffer = ref.load_buffer();

    // 1. Rotate Left
    // Expected: {1, 2, 3, 0} (for Width 4)
    buffer.rotate_left();

    auto pos_x = buffer.position.x.to_array();
    constexpr size_t Width = packedd::size();

    for(size_t i = 0; i < Width; ++i) {
        // The value at index i should be the value originally at (i + 1)
        // Wraps around at the end: index[Width-1] gets value 0.0
        const double expected = static_cast<double>((i + 1) % Width);

        EXPECT_DOUBLE_EQ(pos_x[i], expected)
            << "Rotate Left failed at lane " << i;
    }

    // 2. Rotate Right (Should restore original state)
    // Expected: {0, 1, 2, 3}
    buffer.rotate_right();

    pos_x = buffer.position.x.to_array();
    for(size_t i = 0; i < Width; ++i) {
        const double expected = static_cast<double>(i);
        EXPECT_DOUBLE_EQ(pos_x[i], expected)
            << "Rotate Right (Restore) failed at lane " << i;
    }
}

// 8. Explicit Commit (Write Back)
// Verifies that we can write specific fields back to memory for ALL lanes.
TEST_F(PackedParticleViewsTest, BufferCommit) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask, NoParticleAttributes> ref(src);

    auto buffer = ref.load_buffer();

    // Accumulate force in registers (Memory 0.0 -> Buffer 10.0)
    buffer.force += pvec3(10.0, 10.0, 10.0);

    // Modify position in buffer (Memory i -> Buffer 999.0)
    buffer.position = pvec3(999.0);

    // COMMIT only Force
    ref.force = buffer.force;

    // Verify Memory Updates
    constexpr size_t Width = packedd::size();

    // Check backing vectors directly to be sure
    for(size_t i = 0; i < Width; ++i) {
        // Force should be updated to 10.0 everywhere
        EXPECT_DOUBLE_EQ(force_x[i], 10.0) << "Force commit failed at lane " << i;

        // Position should remain UNTOUCHED (still equal to index i)
        // It must NOT be 999.0
        const double expected_pos = static_cast<double>(i);
        EXPECT_DOUBLE_EQ(pos_x[i], expected_pos) << "Position accidentally overwritten at lane " << i;
    }
}

// 9. Systolic Simulation
// A robust integration test for the loop logic.
TEST_F(PackedParticleViewsTest, SystolicLoopSim) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask, NoParticleAttributes> p1(src);

    // p1 and p2 point to same memory: PosX = {0, 1, 2, 3}
    const PackedParticleRef<TestMask, TestMask, NoParticleAttributes> p2(src);

    auto b1 = p1.load_buffer();
    auto b2 = p2.load_buffer();

    // Step 1: Rotate b2 Left
    // b1 PosX: {0, 1, 2, 3}
    // b2 PosX: {1, 2, 3, 0}
    b2.rotate_left();

    // Interaction: F = (pos2 - pos1)
    const auto delta = b2.position - b1.position;

    // Accumulate
    b1.force += delta;

    // Commit p1 force
    p1.force = b1.force;

    // Verification Logic:
    // Lane 0: 1.0 - 0.0 = 1.0
    // Lane 1: 2.0 - 1.0 = 1.0
    // ...
    // Lane N (Last): 0.0 - N = -N

    constexpr size_t Width = packedd::size();

    for(size_t i = 0; i < Width; ++i) {
        double expected_force;

        if (i < Width - 1) {
            expected_force = 1.0;
        } else {
            // The wrap-around case at the end of the vector
            expected_force = -static_cast<double>(Width - 1);
        }

        EXPECT_DOUBLE_EQ(force_x[i], expected_force)
            << "Systolic math mismatch at lane " << i;
    }
}


// 10. Broadcast from Scalar
// Verifies that a packed buffer is correctly populated across all lanes
// when constructed from a single scalar particle accessor.
TEST_F(PackedParticleViewsTest, BufferBroadcastFromScalar) {
    // 1. Setup a single scalar particle
    double s_pos_x = 1.5, s_pos_y = 2.5, s_pos_z = 3.5;
    double s_vel_x = 4.5, s_vel_y = 5.5, s_vel_z = 6.5;
    double s_frc_x = 7.5, s_frc_y = 8.5, s_frc_z = 9.5;
    double s_mass  = 10.5;

    ParticleSource<TestMask, TestMask, NoParticleAttributes> scalar_src;
    scalar_src.position.x = &s_pos_x;
    scalar_src.position.y = &s_pos_y;
    scalar_src.position.z = &s_pos_z;
    scalar_src.velocity.x = &s_vel_x;
    scalar_src.velocity.y = &s_vel_y;
    scalar_src.velocity.z = &s_vel_z;
    scalar_src.force.x = &s_frc_x;
    scalar_src.force.y = &s_frc_y;
    scalar_src.force.z = &s_frc_z;
    scalar_src.mass = &s_mass;

    // Create the scalar reference
    ScalarParticleRef<TestMask, TestMask, NoParticleAttributes> scalar_ref(scalar_src);

    // 2. Broadcast into a SIMD buffer
    auto buffer = scalar_ref.broadcast();

    // 3. Verify all lanes contain the exact scalar values
    constexpr size_t Width = packedd::size();

    auto pos_x = buffer.position.x.to_array();
    auto pos_y = buffer.position.y.to_array();
    auto pos_z = buffer.position.z.to_array();

    auto vel_x = buffer.velocity.x.to_array();
    auto vel_y = buffer.velocity.y.to_array();
    auto vel_z = buffer.velocity.z.to_array();

    auto frc_x = buffer.force.x.to_array();
    auto frc_y = buffer.force.y.to_array();
    auto frc_z = buffer.force.z.to_array();

    auto mass_arr = buffer.mass.to_array();

    for(size_t i = 0; i < Width; ++i) {
        EXPECT_DOUBLE_EQ(pos_x[i], 1.5) << "Broadcast failed on position.x at lane " << i;
        EXPECT_DOUBLE_EQ(pos_y[i], 2.5) << "Broadcast failed on position.y at lane " << i;
        EXPECT_DOUBLE_EQ(pos_z[i], 3.5) << "Broadcast failed on position.z at lane " << i;

        EXPECT_DOUBLE_EQ(vel_x[i], 4.5) << "Broadcast failed on velocity.x at lane " << i;
        EXPECT_DOUBLE_EQ(vel_y[i], 5.5) << "Broadcast failed on velocity.y at lane " << i;
        EXPECT_DOUBLE_EQ(vel_z[i], 6.5) << "Broadcast failed on velocity.z at lane " << i;

        EXPECT_DOUBLE_EQ(frc_x[i], 7.5) << "Broadcast failed on force.x at lane " << i;
        EXPECT_DOUBLE_EQ(frc_y[i], 8.5) << "Broadcast failed on force.y at lane " << i;
        EXPECT_DOUBLE_EQ(frc_z[i], 9.5) << "Broadcast failed on force.z at lane " << i;

        EXPECT_DOUBLE_EQ(mass_arr[i], 10.5) << "Broadcast failed on mass at lane " << i;
    }
}


// ---------------------
// BUFFER AND VIEW TESTS
// ---------------------

// Unmasked update_into (Write-Back)
// Verifies that modifying a buffer and calling update_into correctly flushes
// registers back to memory across all lanes.
TEST_F(PackedParticleViewsTest, BufferUpdateIntoUnmasked) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask, NoParticleAttributes> ref(src);

    auto buffer = ref.load_buffer();

    // Modify buffer registers
    buffer.position += pvec3(5.0);
    buffer.force = pvec3(100.0);

    // Flush back to memory
    buffer.update_into(ref);

    constexpr size_t Width = packed::size();
    for (size_t i = 0; i < Width; ++i) {
        // Original pos was 'i', should now be 'i + 5.0'
        EXPECT_DOUBLE_EQ(pos_x[i], static_cast<double>(i) + 5.0);
        // Force should be overwritten to 100.0
        EXPECT_DOUBLE_EQ(force_x[i], 100.0);
    }
}

// Masked update_into (Tail Chunk Safety)
// Verifies that update_into respects SIMD masks, preventing memory corruption
// in lanes that are outside the active bounds.
TEST_F(PackedParticleViewsTest, BufferUpdateIntoMasked) {
    const auto src = get_source();
    PackedParticleRef<TestMask, TestMask, NoParticleAttributes> ref(src);

    auto buffer = ref.load_buffer();

    // Modify buffer registers
    buffer.position += pvec3(50.0);

    // Create a mask that is only TRUE for the first 2 lanes
    constexpr size_t Width = packed::size();
    std::vector<double> seq(Width);
    for(size_t i = 0; i < Width; ++i) seq[i] = static_cast<double>(i);

    packed indices = packed::load_unaligned(seq.data());
    auto mask = (indices < 2.0);

    // Masked flush
    buffer.update_into(ref, mask);

    for (size_t i = 0; i < Width; ++i) {
        if (i < 2) {
            // First two lanes should be updated
            EXPECT_DOUBLE_EQ(pos_x[i], static_cast<double>(i) + 50.0);
        } else {
            // Remaining lanes MUST remain untouched
            EXPECT_DOUBLE_EQ(pos_x[i], static_cast<double>(i));
        }
    }
}

// Buffer to View Const-Correctness
// Verifies that to_view() correctly maps read-only fields and exposes them
// without allowing assignment.
TEST_F(PackedParticleViewsTest, BufferToView) {
    const auto src = get_source();
    PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes> ref(src);

    auto buffer = ref.load_buffer();
    auto view = buffer.to_view();

    // 1. Reading is allowed
    const auto v_x = view.position.x;
    EXPECT_DOUBLE_EQ(v_x.to_array()[0], 0.0);

    // 2. Writing must fail to compile
    // We check if we can assign a 'packed' register to the view's position.x
    static_assert(!std::is_assignable_v<decltype((view.position.x)), packed>,
        "FATAL: BufferView is allowing mutable assignments to its fields!");
}

// Masked Reduce Into (Scalar Write-Back)
// Verifies that a pure Write-Only buffer can horizontally reduce its SIMD lanes
// into a single scalar particle (used during Newton's 3rd Law / Cell sorting)
TEST(PackedParticleReductionTest, MaskedReduceIntoScalar) {
    // Setup a pure Write-Only target
    constexpr auto ROMask = ParticleField::none;
    constexpr auto WOMask = ParticleField::force;

    // Single scalar particle to receive the reduced forces
    double s_frc_x = 0.0, s_frc_y = 0.0, s_frc_z = 0.0;
    ParticleSource<ROMask, WOMask, NoParticleAttributes> s_src;
    s_src.force.x = &s_frc_x;
    s_src.force.y = &s_frc_y;
    s_src.force.z = &s_frc_z;
    ScalarParticleRef<ROMask, WOMask, NoParticleAttributes> scalar_ref(s_src);

    // Create a dummy buffer and accumulate some forces into it
    PackedParticleBuffer<ROMask, WOMask, NoParticleAttributes> buffer;
    buffer.force = pvec3(1.0, 2.0, 3.0); // Every lane has {1, 2, 3}

    // Create a mask that is only TRUE for the first 3 lanes
    constexpr size_t Width = packed::size();
    std::vector<double> seq(Width);
    for(size_t i = 0; i < Width; ++i) seq[i] = static_cast<double>(i);

    packed indices = packed::load_unaligned(seq.data());
    auto mask = (indices < 2.0);

    // Reduce
    buffer.reduce_into(scalar_ref, mask);

    // Expected: 2 lanes * 1.0 = 3.0
    EXPECT_DOUBLE_EQ(s_frc_x, 2.0);
    EXPECT_DOUBLE_EQ(s_frc_y, 4.0);
    EXPECT_DOUBLE_EQ(s_frc_z, 6.0);
}


//----------------
// ATTRIBUTE TESTS
//----------------
struct TestCharge {
    APRIL_TRIVIAL_ATTRIBUTE(double, q);
};

// Trivial SIMD Attributes Lifecycle
// Verifies the engine correctly casts and vectorizes custom user structs
TEST(PackedParticleAttributesTest, SIMDAttributeLifecycle) {
    constexpr size_t Count = packed::size();

    // Memory setup
    std::vector<double> pos_x(Count, 0.0);
    std::vector<double> pos_y(Count, 0.0);
    std::vector<double> pos_z(Count, 0.0);
    std::vector<TestCharge> charges(Count);
    for(size_t i=0; i<Count; ++i) {
        charges[i].q = static_cast<double>(i) * 2.0; // 0.0, 2.0, 4.0, ...
    }

    constexpr auto Mask = ParticleField::position | ParticleField::attributes;
    ParticleSource<Mask, Mask, TestCharge> src;
    src.position.x = pos_x.data();
    src.position.y = pos_y.data();
    src.position.z = pos_z.data();

    src.attributes = charges.data();

    // Load
    PackedParticleRef<Mask, Mask, TestCharge> ref(src);
    auto buffer = ref.load_buffer();

    // Verify Read (Via View mapping)
    auto view = buffer.to_view();
    auto q_vals = view.attributes.q.to_array(); // Beautiful user syntax

    for(size_t i=0; i<Count; ++i) {
        EXPECT_DOUBLE_EQ(q_vals[i], static_cast<double>(i) * 2.0);
    }

    // Modify & Write-back
    view.attributes.q += view.attributes.q; // Double the charges
    buffer.update_into(ref);

    for(size_t i=0; i<Count; ++i) {
        EXPECT_DOUBLE_EQ(buffer.attributes.to_array()[i], static_cast<double>(i) * 4.0);
    }

    // Verify Memory
    for(size_t i=0; i<Count; ++i) {
        EXPECT_DOUBLE_EQ(charges[i].q, static_cast<double>(i) * 4.0);
    }
}









