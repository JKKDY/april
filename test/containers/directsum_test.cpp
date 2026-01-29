#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::AnyOf;
using testing::Eq;

#include "april/base/types.hpp"
#include "constant_force.h"
#include "utils.h"

#include "april/containers/direct_sum/ds_aos.hpp"
#include "april/containers/direct_sum/ds_soa.hpp"
#include "april/containers/direct_sum/ds_aosoa.hpp"

using namespace april;

template <typename T>
class DirectSumTest : public testing::Test {};

using ContainerTypes = testing::Types<DirectSumAoS, DirectSumSoA, DirectSumAoSoA>;
TYPED_TEST_SUITE(DirectSumTest, ContainerTypes);

TYPED_TEST(DirectSumTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
	e.add_particle(make_particle(0, {1,2,3}, {}, 1, ParticleState::ALIVE, 0));
	e.add_force(NoForce(), to_type(0));
	e.set_extent(1,1,1);

	auto sys = build_system(e, DirectSumAoS());
    sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TYPED_TEST(DirectSumTest, TwoParticles_ConstantTypeForce) {
    Environment e (forces<ConstantForce>);
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1,0,0}, {}, 1, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));
	e.set_extent(1,1,1);

	auto sys = build_system(e, TypeParam());
    sys.update_forces();
    auto const& out = export_particles(sys);

    ASSERT_EQ(out.size(), 2u);
    // both should see the same force vector
    EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5)))
	);

}

TYPED_TEST(DirectSumTest, TwoParticles_IdSpecificForce) {
    Environment e (forces<ConstantForce, NoForce>);
	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 42));
	e.add_particle(make_particle(0, {0,1,0}, {}, 1, ParticleState::ALIVE, 99));
	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));
	e.set_extent(1,1,1);

	auto sys = build_system(e, TypeParam());
	sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TYPED_TEST(DirectSumTest, TwoParticles_InverseSquare) {
	Environment e (forces<Gravity, NoForce>);

	e.set_extent({10,10,10});

	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(1, {2,0,0}, {}, 2, ParticleState::ALIVE, 1));


	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));
	e.add_force(Gravity(5.0), between_types(0, 1));

	auto sys = build_system(e, TypeParam());
    sys.update_forces();

    auto const& out = export_particles(sys);
    // find each
    const auto& pa = (out[0].mass == 1 ? out[0] : out[1]);
    const auto& pb = (out[1].mass == 2 ? out[1] : out[0]);
    // magnitude = pre * m1*m2 / r^3 = 5*1*2/(2^3)=10/8=1.25  direction from pa->pb = (2,0,0)
    // force on pa = 1.25*(2,0,0) = (2.5,0,0); on pb = (-2.5,0,0)
    EXPECT_NEAR(pa.force.x, 2.5, 1e-12);
    EXPECT_NEAR(pb.force.x, -2.5, 1e-12);
    EXPECT_EQ(pa.force.y, 0.0);
    EXPECT_EQ(pb.force.y, 0.0);
}


TYPED_TEST(DirectSumTest, CollectIndicesInRegion) {
	// Create a simple 3x3x3 grid of particles (27 total)
	auto cuboid = ParticleCuboid{}
		.at(vec3(0.25))
		.velocity({0, 0, 0})
		.count({3, 3, 3})
		.mass(1.0)
		.spacing(1)
		.type(0);

    Environment e(forces<NoForce>);
    e.set_origin({0, 0, 0});
    e.set_extent({5, 5, 5});
	e.add_particles(cuboid);
    e.add_force(NoForce(), to_type(0));

    auto sys = build_system(e, TypeParam());

    // Case 1: small inner region (should include one particle)
    {
        env::Domain region({0.1, 0.1, 0.1}, {0.9, 0.9, 0.9});
        auto indices = sys.query_region(env::Box::from_domain(region));
        ASSERT_EQ(indices.size(), 1u);
        auto p = export_particles(sys)[indices[0]];
        EXPECT_EQ(p.position.x, 0.25);
        EXPECT_EQ(p.position.y, 0.25);
        EXPECT_EQ(p.position.z, 0.25);
    }

    // Case 2: mid region (should include all 27)
    {
        env::Domain region({0, 0, 0}, {5, 5, 5});
        auto indices = sys.query_region(env::Box::from_domain(region));
        EXPECT_EQ(indices.size(), 27u);
    }

    // Case 3: partially overlapping region
    {
        env::Domain region({1.5, 1.5, 1.5}, {4.5, 4.5, 4.5});
        std::vector indices = sys.query_region(env::Box::from_domain(region));
        EXPECT_GT(indices.size(), 0u);
        EXPECT_LT(indices.size(), 27u);

        std::unordered_set inside (indices.begin(), indices.end());

        for (size_t id = sys.min_id(); id < sys.max_id(); ++id) {
			auto p = get_particle(sys, id);
        	bool in_region = (p.position.x >= 1.5 && p.position.x <= 4.5) &&
				(p.position.y >= 1.5 && p.position.y <= 4.5) &&
				(p.position.z >= 1.5 && p.position.z <= 4.5);

        	if (inside.contains(id)) {
        		EXPECT_TRUE(in_region);
        	} else {
        		EXPECT_FALSE(in_region);
        	}
        }
    }

    // Case 4: region completely outside
    {
        env::Domain region({10, 10, 10}, {12, 12, 12});
        auto indices = sys.query_region(env::Box::from_domain(region));
        EXPECT_TRUE(indices.empty());
    }
}


// does nothing except signaling the container to be periodic
struct DummyPeriodicBoundary final : Boundary {
	static constexpr env::FieldMask fields = +env::Field::none;

	DummyPeriodicBoundary()
	: Boundary(0.0, false, true, false ) {}

	template<env::FieldMask M, env::IsUserData U>
	void apply(env::ParticleRef<M, U> &, const env::Box &, const Face) const noexcept {}
};

TYPED_TEST(DirectSumTest, PeriodicForceWrap_X) {
	Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);

	e.set_origin({0,0,0});
	e.set_extent({10,10,10}); // domain box 10x10x10

	// Two particles, near opposite faces along x
	e.add_particle(make_particle(0, {0.5, 5, 5}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(0, {9.5, 5, 5}, {}, 1, ParticleState::ALIVE, 1));

	e.add_force(Harmonic(1, 0, 2), to_type(0)); // simple directional force
	e.set_boundaries(DummyPeriodicBoundary(), {Face::XMinus, Face::XPlus});

	BuildInfo mapping;
	auto sys = build_system(e, TypeParam(),&mapping); // DirectSum container
	sys.update_forces();

	const auto out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
	auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

	// They should feel equal and opposite forces due to wrapping
	EXPECT_EQ(p1.force, -p2.force);
	EXPECT_EQ(p1.force.x, -1.0);
	EXPECT_EQ(p2.force.x, 1.0);
}


TYPED_TEST(DirectSumTest, PeriodicForceWrap_AllAxes) {
	// Enable force wrapping in all 6 directions
	Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
	e.set_origin({0, 0, 0});
	e.set_extent({10, 10, 10});

	// Particles placed in diagonally opposite corners
	e.add_particle(make_particle(0, {0.5, 0.5, 0.5}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(0, {9.5, 9.5, 9.5}, {}, 1, ParticleState::ALIVE, 1));

	// Hooke-like spring with k=1, r0=0, cutoff=2
	e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

	// Activate periodic wrapping on all faces
	e.set_boundaries(DummyPeriodicBoundary(), {
		Face::XMinus, Face::XPlus,
		Face::YMinus, Face::YPlus,
		Face::ZMinus, Face::ZPlus
	});

	BuildInfo mapping;
	auto sys = build_system(e, TypeParam(), &mapping);
	sys.update_forces();

	auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
	auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

	// In a 10x10x10 domain with full wrapping,
	// the wrapped displacement should be (-1, -1, -1) for p1->p2.
	// The harmonic force acts along that direction, magnitude proportional to distance.

	// Check that forces are opposite
	EXPECT_EQ(p1.force, -p2.force);

	// Check that the direction is consistent with wrapped displacement
	EXPECT_EQ(p1.force.x, -1.0);
	EXPECT_EQ(p1.force.y, -1.0);
	EXPECT_EQ(p1.force.z, -1.0);

	EXPECT_EQ(p2.force.x, 1.0);
	EXPECT_EQ(p2.force.y, 1.0);
	EXPECT_EQ(p2.force.z, 1.0);
}


TYPED_TEST(DirectSumTest, Asymmetric_ChunkBoundaries_Counting) {
    constexpr size_t n_type0 = 20;
    constexpr size_t n_type1 = 12;

    Environment e(forces<ConstantForce, NoForce>); // Added NoForce
    e.set_extent({10, 10, 10});

    for (ParticleID i = 0; i < n_type0; ++i) {
        e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, i));
    }
    for (ParticleID i = 0; i < n_type1; ++i) {
        e.add_particle(make_particle(1, {1,0,0}, {}, 1, ParticleState::ALIVE, 100 + i));
    }

    e.add_force(ConstantForce(1, 2, 3), between_types(0, 1));

    // Explicitly define self-interactions as NoForce
    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));

	BuildInfo info;
    auto sys = build_system(e, TypeParam(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), n_type0 + n_type1);

    const vec3 expected_f0 = vec3(-1, -2, -3) * static_cast<double>(n_type1);
    const vec3 expected_f1 = vec3(1, 2, 3) * static_cast<double>(n_type0);

	for (const auto& p : out) {
		if (p.type == info.type_map[0]) {
			EXPECT_THAT(p.force, testing::AnyOf(expected_f0, -expected_f0));
		} else {
			EXPECT_THAT(p.force, testing::AnyOf(expected_f1, -expected_f1));
		}
	}
}


TYPED_TEST(DirectSumTest, Asymmetric_MultiChunk_Gravity) {
    constexpr size_t n_a = 10;
    constexpr size_t n_b = 10;

    Environment e(forces<Gravity, NoForce>);
    e.set_extent({100, 100, 100});

    for (ParticleID i = 0; i < n_a; ++i) {
        e.add_particle(make_particle(0, {static_cast<double>(i), 0, 0}, {}, 1.0, ParticleState::ALIVE, i));
    }
    for (ParticleID i = 0; i < n_b; ++i) {
        e.add_particle(make_particle(1, {static_cast<double>(i), 10, 0}, {}, 1.0, ParticleState::ALIVE, 100+i));
    }

    e.add_force(Gravity(1.0), between_types(0, 1));
    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));

    BuildInfo info;
    auto sys = build_system(e, TypeParam(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);

    // Verify against O(N^2) brute force
    for (const auto& p : out) {
        vec3 expected_force(0,0,0);
        // Map internal type back to user type for logic check, or compare against mapped type
        const int target_user_type = (p.type == info.type_map[0]) ? 1 : 0;

        for (const auto& other : out) {
            if (other.type != info.type_map[target_user_type]) continue;

            vec3 r = other.position - p.position;
            const double dist = r.norm();
            const double mag = (1.0 * p.mass * other.mass) / (dist * dist * dist);
            expected_force += r * mag;
        }

        EXPECT_NEAR(p.force.x, expected_force.x, 1e-10);
        EXPECT_NEAR(p.force.y, expected_force.y, 1e-10);
        EXPECT_NEAR(p.force.z, expected_force.z, 1e-10);
    }
}

TYPED_TEST(DirectSumTest, Asymmetric_TypeChaining) {
    Environment e(forces<Harmonic, NoForce>);
    e.set_extent({10, 10, 10});

    // P0(0,0,0) --[k=100]--> P1(1,0,0) --[k=10]--> P2(1,1,0)
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
    e.add_particle(make_particle(1, {1,0,0}, {}, 1, ParticleState::ALIVE, 1));
    e.add_particle(make_particle(2, {1,1,0}, {}, 1, ParticleState::ALIVE, 2));

    e.add_force(Harmonic(100, 0, 5), between_types(0, 1));
    e.add_force(Harmonic(10, 0, 5),  between_types(1, 2));

    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));
    e.add_force(NoForce(), to_type(2));
    e.add_force(NoForce(), between_types(0, 2));

    BuildInfo info;
    auto sys = build_system(e, TypeParam(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);

    auto get_force = [&](const size_t id) {
        for (const auto& p : out) if (p.id == info.id_map[id]) return p.force;
        return vec3(0);
    };

    // P0 is pulled +X towards P1
    EXPECT_EQ(get_force(0), vec3(100, 0, 0));

    // P2 is pulled -Y towards P1
    EXPECT_EQ(get_force(2), vec3(0, -10, 0));

    // P1 is pulled -X towards P0 and +Y towards P2
    EXPECT_EQ(get_force(1), vec3(-100, 10, 0));
}