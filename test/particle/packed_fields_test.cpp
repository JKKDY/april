#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>
#include <concepts>


#include "april/particle/access/source.hpp"
#include "april/particle/access/packed_access.hpp"
#include "april/base/types.hpp"

using namespace april;
using namespace april::particle;
using namespace april::particle::internal;

// Define a test mask (Standard Physics Fields)
static constexpr auto TestMask = ParticleField::position | ParticleField::velocity | ParticleField::force | ParticleField::mass;

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

        for (size_t i = 0; i < Count; ++i) {
            const double val = static_cast<double>(i);

            pos_x[i] = val;
            pos_y[i] = val;
            pos_z[i] = val;

            vel_x[i] = 1.0;
            vel_y[i] = 0.0;
            vel_z[i] = 0.0;

            mass[i] = 2.0;
        }
    }

public:
    [[nodiscard]] auto get_source() {
        auto get_field = [&]<ParticleField F>() {
            if constexpr (F == ParticleField::position)
                return math::Vec3Location{pos_x.data(), pos_y.data(), pos_z.data()};
            else if constexpr (F == ParticleField::velocity)
                return math::Vec3Location{vel_x.data(), vel_y.data(), vel_z.data()};
            else if constexpr (F == ParticleField::force)
                return math::Vec3Location{force_x.data(), force_y.data(), force_z.data()};
            else if constexpr (F == ParticleField::mass)
                return mass.data();
        };

        return make_particle_source<TestMask, TestMask>(get_field);
    }

    void ExpectParticle(
        const size_t index,
        const vec3& expected_pos,
        const vec3& expected_force
    ) const {
        EXPECT_DOUBLE_EQ(pos_x[index], expected_pos.x);
        EXPECT_DOUBLE_EQ(pos_y[index], expected_pos.y);
        EXPECT_DOUBLE_EQ(pos_z[index], expected_pos.z);

        EXPECT_DOUBLE_EQ(force_x[index], expected_force.x);
        EXPECT_DOUBLE_EQ(force_y[index], expected_force.y);
        EXPECT_DOUBLE_EQ(force_z[index], expected_force.z);
    }
};

using Source = decltype(std::declval<PackedParticleViewsTest&>().get_source());
using PackedParticle = PackedParticleRef<TestMask, TestMask, NoParticleAttributes, Source>;

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
    PackedParticleRef<TestMask, TestMask,NoParticleAttributes, decltype(src)> ref(src);

    // Convert to View
    const PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes, decltype(src)> view = ref.to_view();

    // Read is allowed
    const pvec3 v = view.velocity;
    EXPECT_DOUBLE_EQ(v.x.to_array()[0], 1.0);

    // Write should fail compile on view, should succeed on ref
    // TODO
    // static_assert(!VelocityAssignable<PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes>>);
    static_assert( VelocityAssignable<PackedParticleRef<TestMask, TestMask, NoParticleAttributes, decltype(src)>>);
}



// 6. Buffer Load & Isolation
// Verifies that 'load_buffer' reads correctly across ALL lanes,
// and that modifying the registers DOES NOT touch memory until commit.
TEST_F(PackedParticleViewsTest, BufferIsolation) {
    const auto src = get_source();
    PackedParticle ref(src);

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
    const PackedParticle ref(src);

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
    PackedParticle ref(src);

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
    PackedParticle p1(src);

    // p1 and p2 point to same memory: PosX = {0, 1, 2, 3}
    const PackedParticle p2(src);

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
	double s_mass = 10.5;

	auto get_field = [&]<ParticleField F>() {
		if constexpr (F == ParticleField::position)
			return math::Vec3Location{&s_pos_x, &s_pos_y, &s_pos_z};
		else if constexpr (F == ParticleField::velocity)
			return math::Vec3Location{&s_vel_x, &s_vel_y, &s_vel_z};
		else if constexpr (F == ParticleField::force)
			return math::Vec3Location{&s_frc_x, &s_frc_y, &s_frc_z};
		else if constexpr (F == ParticleField::mass)
			return &s_mass;
	};

	auto scalar_src = particle::internal::make_particle_source<TestMask, TestMask>(get_field);
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

	for (size_t i = 0; i < Width; ++i) {
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
    PackedParticle ref(src);

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
    PackedParticle ref(src);

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
    PackedParticleRef<TestMask, ParticleField::none, NoParticleAttributes, decltype(src)> ref(src);

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
    constexpr auto ROMask = ParticleField::none;
    constexpr auto WOMask = ParticleField::force;

    double s_frc_x = 0.0, s_frc_y = 0.0, s_frc_z = 0.0;

    auto get_field = [&]<ParticleField F>() {
        return math::Vec3Location{
            &s_frc_x,
            &s_frc_y,
            &s_frc_z
        };
    };

    auto s_src = make_particle_source<ROMask, WOMask>(get_field);
    ScalarParticleRef<ROMask, WOMask, NoParticleAttributes> scalar_ref(s_src);

    PackedParticleBuffer<
        ROMask,
        WOMask,
        NoParticleAttributes,
        MaskPolicy::Disabled
    > buffer;

    buffer.force = pvec3(1.0, 2.0, 3.0);

    constexpr size_t Width = packed::size();
    std::vector<double> seq(Width);
    for (size_t i = 0; i < Width; ++i) seq[i] = static_cast<double>(i);

    packed indices = packed::load_unaligned(seq.data());
    auto mask = indices < 2.0;

    buffer.reduce_into(scalar_ref, mask);

    EXPECT_DOUBLE_EQ(s_frc_x, 2.0);
    EXPECT_DOUBLE_EQ(s_frc_y, 4.0);
    EXPECT_DOUBLE_EQ(s_frc_z, 6.0);
}


//----------------
// ATTRIBUTE TESTS
//----------------





struct TestAttributes {
    double charge = 0.0;
    size_t touched_by_lane = std::numeric_limits<size_t>::max();
};

TEST(PackedParticleAttributesTest, AttributesAreForwardedIntoPackedKernel) {
    constexpr size_t Width = packed::size();

    std::array<TestAttributes, Width> attributes{};

    for (size_t lane = 0; lane < Width; ++lane)
        attributes[lane].charge = static_cast<double>(lane + 1);

    constexpr ParticleField Read =
        ParticleField::attributes;

    constexpr ParticleField Write =
        ParticleField::attributes;

    auto get_field = [&]<ParticleField F>() {
        if constexpr (F == ParticleField::attributes)
            return attributes.data();
    };

    auto source = make_particle_source<Read, Write>(get_field);

    using AttributeSource = decltype(source);

    using Ref = PackedParticleRef<
        Read,
        Write,
        TestAttributes,
        AttributeSource
    >;

    Ref ref(source);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    std::array<bool, Width> active_lanes{};

    for (size_t lane = 0; lane < Width; ++lane)
        active_lanes[lane] = (lane % 2) == 0;

    const auto mask =
        packed::mask_type::load_unaligned(active_lanes.data());

    view.mask_with(mask);

    // Simulated packed kernel
    for (size_t lane = 0; lane < Width; ++lane) {
        if (!view.attributes.active(lane))
            continue;

        auto& attr = view.attributes[lane];

        attr.charge *= 10.0;
        attr.touched_by_lane = lane;
    }

    // Attribute writes are direct; no update_into().
    for (size_t lane = 0; lane < Width; ++lane) {
        if ((lane % 2) == 0) {
            EXPECT_DOUBLE_EQ(
                attributes[lane].charge,
                static_cast<double>(lane + 1) * 10.0
            );

            EXPECT_EQ(
                attributes[lane].touched_by_lane,
                lane
            );
        } else {
            EXPECT_DOUBLE_EQ(
                attributes[lane].charge,
                static_cast<double>(lane + 1)
            );

            EXPECT_EQ(
                attributes[lane].touched_by_lane,
                std::numeric_limits<size_t>::max()
            );
        }
    }
}

TEST(PackedParticleAttributesTest, ReadOnlyAttributesAreConstInKernel) {
    constexpr size_t Width = packed::size();

    std::array<TestAttributes, Width> attributes{};

    constexpr ParticleField Read =
        ParticleField::attributes;

    constexpr ParticleField Write =
        ParticleField::none;

    auto get_field = [&]<ParticleField F>() {
        if constexpr (F == ParticleField::attributes)
            return attributes.data();
    };

    auto source = make_particle_source<Read, Write>(get_field);

    using AttributeSource = decltype(source);

    using Ref = PackedParticleRef<
        Read,
        Write,
        TestAttributes,
        AttributeSource
    >;

    Ref ref(source);

    auto buffer = ref.load_buffer();
    auto view = buffer.to_view();

    using AttributeRef =
        decltype(view.attributes[0]);

    static_assert(std::same_as<
        AttributeRef,
        const TestAttributes&
    >);

    EXPECT_EQ(
        &view.attributes[0],
        &attributes[0]
    );
}


namespace {

using ParticleMaskT = packed::mask_type;

template<typename Predicate>
ParticleMaskT MakeParticleMask(Predicate&& predicate) {
    std::array<bool, packed::size()> lanes{};

    for (std::size_t i = 0; i < lanes.size(); ++i) {
        lanes[i] = static_cast<bool>(predicate(i));
    }

    return ParticleMaskT::load_unaligned(lanes.data());
}

ParticleMaskT FirstParticleLanes(std::size_t count) {
    return MakeParticleMask([count](std::size_t lane) {
        return lane < count;
    });
}

ParticleMaskT AlternatingParticleLanes(
    bool first_lane_active = true
) {
    return MakeParticleMask([first_lane_active](std::size_t lane) {
        return ((lane % 2) == 0) == first_lane_active;
    });
}

template<typename View>
concept SupportsMaskWith = requires(
    View& view,
    const ParticleMaskT& mask
) {
    view.mask_with(mask);
};

} // namespace


// ---------------------------------------------------------
// COMPILE-TIME MASKING POLICY
// ---------------------------------------------------------


TEST(PackedParticleMaskedBufferTest, OnlyWritableFieldsAreMasked) {
    constexpr ParticleField Read =
        ParticleField::position |
        ParticleField::mass;

    constexpr ParticleField Write =
        ParticleField::force;

    using Buffer = PackedParticleBuffer<
        Read,
        Write,
        NoParticleAttributes,
        MaskPolicy::Enabled
    >;

    using PositionComponent = std::remove_cvref_t<
        decltype(std::declval<Buffer&>().position.x)
    >;

    using ForceComponent = std::remove_cvref_t<
        decltype(std::declval<Buffer&>().force.x)
    >;

    using MassField = std::remove_cvref_t<
        decltype(std::declval<Buffer&>().mass)
    >;

    static_assert(std::same_as<PositionComponent, packed>);
    static_assert(std::same_as<MassField, packed>);

    static_assert(std::same_as<
        ForceComponent,
        simd::MaskedPacked<packed, ParticleMaskT>
    >);
}


TEST(PackedParticleMaskedBufferTest, MaskWithExistsOnlyOnMaskedView) {
    using MaskedBuffer = PackedParticleBuffer<
        TestMask,
        TestMask,
        NoParticleAttributes,
        MaskPolicy::Enabled
    >;

    using UnmaskedBuffer = PackedParticleBuffer<
        TestMask,
        TestMask,
        NoParticleAttributes,
        MaskPolicy::Disabled
    >;

    using MaskedView = decltype(
        std::declval<MaskedBuffer&>().to_view()
    );

    using UnmaskedView = decltype(
        std::declval<UnmaskedBuffer&>().to_view()
    );

    static_assert(SupportsMaskWith<MaskedView>);
    static_assert(!SupportsMaskWith<UnmaskedView>);
}


// ---------------------------------------------------------
// DEFAULT MASK STATE
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    EnabledBufferInitiallyAllowsWritesToEveryLane
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    view.position.x += 5.0;

    buffer.update_into(ref);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        EXPECT_DOUBLE_EQ(
            pos_x[i],
            static_cast<double>(i) + 5.0
        ) << "lane " << i;
    }
}


// ---------------------------------------------------------
// READ-WRITE FIELD MASKING
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    MaskedViewUpdatesOnlyActiveReadWriteLanes
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    const auto active = AlternatingParticleLanes();
    view.mask_with(active);

    view.position += pvec3(10.0, 20.0, 30.0);
    view.velocity = pvec3(7.0, 8.0, 9.0);
    view.force = pvec3(4.0, 5.0, 6.0);
    view.mass *= 3.0;

    buffer.update_into(ref);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        const bool lane_active = (i % 2) == 0;
        const double initial_position = static_cast<double>(i);

        EXPECT_DOUBLE_EQ(
            pos_x[i],
            lane_active
                ? initial_position + 10.0
                : initial_position
        ) << "position.x lane " << i;

        EXPECT_DOUBLE_EQ(
            pos_y[i],
            lane_active
                ? initial_position + 20.0
                : initial_position
        ) << "position.y lane " << i;

        EXPECT_DOUBLE_EQ(
            pos_z[i],
            lane_active
                ? initial_position + 30.0
                : initial_position
        ) << "position.z lane " << i;

        EXPECT_DOUBLE_EQ(
            vel_x[i],
            lane_active ? 7.0 : 1.0
        ) << "velocity.x lane " << i;

        EXPECT_DOUBLE_EQ(
            vel_y[i],
            lane_active ? 8.0 : 0.0
        ) << "velocity.y lane " << i;

        EXPECT_DOUBLE_EQ(
            vel_z[i],
            lane_active ? 9.0 : 0.0
        ) << "velocity.z lane " << i;

        EXPECT_DOUBLE_EQ(
            force_x[i],
            lane_active ? 4.0 : 0.0
        ) << "force.x lane " << i;

        EXPECT_DOUBLE_EQ(
            force_y[i],
            lane_active ? 5.0 : 0.0
        ) << "force.y lane " << i;

        EXPECT_DOUBLE_EQ(
            force_z[i],
            lane_active ? 6.0 : 0.0
        ) << "force.z lane " << i;

        EXPECT_DOUBLE_EQ(
            mass[i],
            lane_active ? 6.0 : 2.0
        ) << "mass lane " << i;
    }
}


// ---------------------------------------------------------
// PROSPECTIVE AND CUMULATIVE MASKING
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    MaskWithIsProspectiveAndCumulative
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    // The initial write mask is all true.
    view.position.x += 1.0;

    const std::size_t half =
        (packed::size() + 1) / 2;

    // Restrict later mutations to the first half.
    view.mask_with(FirstParticleLanes(half));
    view.position.x += 10.0;

    // mask_with() intersects with the current mask.
    view.mask_with(AlternatingParticleLanes());
    view.position.x += 100.0;

    buffer.update_into(ref);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        double expected = static_cast<double>(i);

        // This mutation occurred before any mask restriction.
        expected += 1.0;

        if (i < half) {
            expected += 10.0;
        }

        if (i < half && (i % 2) == 0) {
            expected += 100.0;
        }

        EXPECT_DOUBLE_EQ(pos_x[i], expected)
            << "lane " << i;
    }
}


// ---------------------------------------------------------
// WRITE-ONLY ACCUMULATION
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    MaskedWriteOnlyFieldAccumulatesOnlyActiveLanes
) {
    constexpr ParticleField Read =
        ParticleField::position |
        ParticleField::velocity |
        ParticleField::mass;

    constexpr ParticleField Write =
        ParticleField::force;

    const auto src = get_source();
    using Ref = PackedParticleRef<
        Read,
        Write,
        NoParticleAttributes,
        Source
    >;

    std::fill(force_x.begin(), force_x.end(), 10.0);
    std::fill(force_y.begin(), force_y.end(), 20.0);
    std::fill(force_z.begin(), force_z.end(), 30.0);

    Ref ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    const auto active = AlternatingParticleLanes();
    view.mask_with(active);

    // A write-only force field begins as a zero delta buffer.
    view.force += pvec3(1.0, 2.0, 3.0);
    view.force += pvec3(4.0, 5.0, 6.0);

    buffer.update_into(ref);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        const bool lane_active = (i % 2) == 0;

        EXPECT_DOUBLE_EQ(
            force_x[i],
            lane_active ? 15.0 : 10.0
        ) << "force.x lane " << i;

        EXPECT_DOUBLE_EQ(
            force_y[i],
            lane_active ? 27.0 : 20.0
        ) << "force.y lane " << i;

        EXPECT_DOUBLE_EQ(
            force_z[i],
            lane_active ? 39.0 : 30.0
        ) << "force.z lane " << i;
    }
}


// ---------------------------------------------------------
// MASK ROTATION
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    BufferRotationKeepsMaskAlignedWithParticleData
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    // Only original lane zero may be mutated.
    view.mask_with(
        MakeParticleMask([](std::size_t lane) {
            return lane == 0;
        })
    );

    // Original lane zero moves to lane one.
    buffer.rotate_right();

    view.position.x += 100.0;

    const auto values =
        static_cast<packed>(view.position.x).to_array();

    constexpr std::size_t Width = packed::size();
    constexpr std::size_t ActiveLane = 1 % Width;

    for (std::size_t i = 0; i < Width; ++i) {
        const std::size_t source_lane =
            (i + Width - 1) % Width;

        double expected =
            static_cast<double>(source_lane);

        if (i == ActiveLane) {
            expected += 100.0;
        }

        EXPECT_DOUBLE_EQ(values[i], expected)
            << "lane " << i;
    }
}


// ---------------------------------------------------------
// INTERNAL AND EXTERNAL MASK COMPOSITION
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    CommitMaskComposesWithBufferWriteMask
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    const std::size_t half =
        (packed::size() + 1) / 2;

    // Controls which lanes are mutated in the register buffer.
    view.mask_with(FirstParticleLanes(half));
    view.position.x += 50.0;

    // Controls which lanes are committed to memory.
    const auto commit_mask =
        AlternatingParticleLanes();

    buffer.update_into(ref, commit_mask);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        const bool changed =
            i < half &&
            (i % 2) == 0;

        const double expected =
            static_cast<double>(i) +
            (changed ? 50.0 : 0.0);

        EXPECT_DOUBLE_EQ(pos_x[i], expected)
            << "lane " << i;
    }
}


// ---------------------------------------------------------
// SHARED MASK ACROSS WRITABLE FIELDS
// ---------------------------------------------------------

TEST_F(
    PackedParticleViewsTest,
    OneViewMaskControlsEveryWritableField
) {
    const auto src = get_source();
    PackedParticle ref(src);

    auto buffer =
        ref.load_buffer<MaskPolicy::Enabled>();

    auto view = buffer.to_view();

    const auto active =
        MakeParticleMask([](std::size_t lane) {
            return lane == 1 || lane == 3;
        });

    view.mask_with(active);

    view.position.x += 100.0;
    view.velocity.y += 200.0;
    view.force.z += 300.0;
    view.mass += 400.0;

    buffer.update_into(ref);

    constexpr std::size_t Width = packed::size();

    for (std::size_t i = 0; i < Width; ++i) {
        const bool lane_active =
            i == 1 || i == 3;

        EXPECT_DOUBLE_EQ(
            pos_x[i],
            static_cast<double>(i) +
                (lane_active ? 100.0 : 0.0)
        ) << "position.x lane " << i;

        EXPECT_DOUBLE_EQ(
            vel_y[i],
            lane_active ? 200.0 : 0.0
        ) << "velocity.y lane " << i;

        EXPECT_DOUBLE_EQ(
            force_z[i],
            lane_active ? 300.0 : 0.0
        ) << "force.z lane " << i;

        EXPECT_DOUBLE_EQ(
            mass[i],
            lane_active ? 402.0 : 2.0
        ) << "mass lane " << i;
    }
}





