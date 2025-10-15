#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::AnyOf;
using testing::Eq;

#include "april/april.h"
using namespace april;

// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final {
	vec3 v;
	double cutoff_radius;
	ConstantForce(double x, double y, double z, double cutoff = -1) : v{x,y,z} {
		cutoff_radius = cutoff;
	}
	vec3 operator()(const env::internal::Particle&, const env::internal::Particle&, const vec3&) const noexcept {
		return v;
	}

	[[nodiscard]] ConstantForce mix(const ConstantForce& other) const noexcept {
		return {
			v.x + other.v.x,
			v.y + other.v.y,
			v.z + other.v.z,
			std::max(cutoff_radius, other.cutoff_radius)
		};
	}
};


TEST(DirectSumTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
    e.add_particle(Particle{.id = 0, .type = 0, .position={1,2,3},.velocity={0,0,0}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));

	auto sys = build_system(e, DirectSum());
    sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(DirectSumTest, TwoParticles_ConstantTypeForce) {
    Environment e (forces<ConstantForce>);
	e.add_particle(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 1, .type = 7, .position={1,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, DirectSum());
    sys.update_forces();
    auto const& out = sys.export_particles();

    ASSERT_EQ(out.size(), 2u);
    // both should see the same force vector
    EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5)))
	);

}

TEST(DirectSumTest, TwoParticles_IdSpecificForce) {
    Environment e (forces<ConstantForce, NoForce>);
    e.add_particle(Particle{.id = 42, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 99, .type = 0, .position={0,1,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));

	auto sys = build_system(e, DirectSum());
	sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TEST(DirectSumTest, TwoParticles_InverseSquare) {
	Environment e (forces<InverseSquare, NoForce>);

	e.set_extent({10,10,10});
	e.add_particle(Particle{.id = 0, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 1, .type = 1, .position={2,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});

	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));
	e.add_force(InverseSquare(5.0), between_types(0, 1));

	auto sys = build_system(e, DirectSum());
    sys.update_forces();

    auto const& out = sys.export_particles();
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


TEST(DirectSumTest, CollectIndicesInRegion) {
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

    auto sys = build_system(e, DirectSum());

    // Case 1: small inner region (should include one particle)
    {
        env::Domain region({0.1, 0.1, 0.1}, {0.9, 0.9, 0.9});
        auto indices = sys.collect_indices_in_region(region);
        ASSERT_EQ(indices.size(), 1u);
        auto p = sys.export_particles()[indices[0]];
        EXPECT_EQ(p.position.x, 0.25);
        EXPECT_EQ(p.position.y, 0.25);
        EXPECT_EQ(p.position.z, 0.25);
    }

    // Case 2: mid region (should include all 27)
    {
        env::Domain region({0, 0, 0}, {5, 5, 5});
        auto indices = sys.collect_indices_in_region(region);
        EXPECT_EQ(indices.size(), 27u);
    }

    // Case 3: partially overlapping region
    {
        env::Domain region({1.5, 1.5, 1.5}, {4.5, 4.5, 4.5});
        std::vector indices = sys.collect_indices_in_region(region);
        EXPECT_GT(indices.size(), 0u);
        EXPECT_LT(indices.size(), 27u);

        std::unordered_set inside (indices.begin(), indices.end());

        for (size_t idx = sys.index_start(); idx < sys.index_end(); idx++) {
			env::internal::Particle & p = sys.get_particle_by_index(idx);
        	bool in_region = (p.position.x >= 1.5 && p.position.x <= 4.5) &&
				(p.position.y >= 1.5 && p.position.y <= 4.5) &&
				(p.position.z >= 1.5 && p.position.z <= 4.5);

        	if (inside.contains(idx)) {
        		EXPECT_TRUE(in_region);
        	} else {
        		EXPECT_FALSE(in_region);
        	}
        }
    }

    // Case 4: region completely outside
    {
        env::Domain region({10, 10, 10}, {12, 12, 12});
        auto indices = sys.collect_indices_in_region(region);
        EXPECT_TRUE(indices.empty());
    }
}


// does nothing except signaling the container to be periodic
struct DummyPeriodicBoundary final : Boundary {
	DummyPeriodicBoundary()
	: Boundary(0.0, false, true, false ) {}

	void apply(env::internal::Particle&, const env::Box&, Face) const noexcept {}
};

TEST(DirectSumTest, PeriodicForceWrap_X) {
	Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);

	e.set_origin({0,0,0});
	e.set_extent({10,10,10}); // domain box 10x10x10

	// Two particles, near opposite faces along x
	e.add_particle(Particle{.id=0, .type=0, .position={0.5, 5, 5}, .velocity={}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_particle(Particle{.id=1, .type=0, .position={9.5, 5, 5}, .velocity={}, .mass=1.0, .state=ParticleState::ALIVE});

	e.add_force(Harmonic(1, 0, 2), to_type(0)); // simple directional force
	e.set_boundaries(DummyPeriodicBoundary(), {Face::XMinus, Face::XPlus});

	UserToInternalMappings mapping;
	auto sys = build_system(e, DirectSum(),&mapping); // DirectSum container
	sys.update_forces();

	auto const& out = sys.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto p1 = sys.get_particle_by_id(mapping.usr_ids_to_impl_ids[0]);
	auto p2 = sys.get_particle_by_id(mapping.usr_ids_to_impl_ids[1]);

	// They should feel equal and opposite forces due to wrapping
	EXPECT_EQ(p1.force, -p2.force);
	EXPECT_EQ(p1.force.x, -1.0);
	EXPECT_EQ(p2.force.x, 1.0);
}


TEST(DirectSumTest, PeriodicForceWrap_AllAxes) {
	// Enable force wrapping in all 6 directions
	Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
	e.set_origin({0, 0, 0});
	e.set_extent({10, 10, 10});

	// Particles placed in diagonally opposite corners
	e.add_particle(Particle{.id=0, .type=0, .position={0.5, 0.5, 0.5}, .velocity={}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_particle(Particle{.id=1, .type=0, .position={9.5, 9.5, 9.5}, .velocity={}, .mass=1.0, .state=ParticleState::ALIVE});

	// Hooke-like spring with k=1, r0=0, cutoff=2
	e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

	// Activate periodic wrapping on all faces
	e.set_boundaries(DummyPeriodicBoundary(), {
		Face::XMinus, Face::XPlus,
		Face::YMinus, Face::YPlus,
		Face::ZMinus, Face::ZPlus
	});

	UserToInternalMappings mapping;
	auto sys = build_system(e, DirectSum(), &mapping);
	sys.update_forces();

	auto const& out = sys.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto p1 = sys.get_particle_by_id(mapping.usr_ids_to_impl_ids[0]);
	auto p2 = sys.get_particle_by_id(mapping.usr_ids_to_impl_ids[1]);

	// In a 10x10x10 domain with full wrapping,
	// the wrapped displacement should be (-1, -1, -1) for p1->p2.
	// The harmonic force acts along that direction, magnitude proportional to distance.

	// Check that forces are opposite
	EXPECT_EQ(p1.force, -p2.force);

	// Check that the direction is consistent with wrapped displacement
	EXPECT_EQ(p1.force.x, -1.0);
	EXPECT_EQ(p1.force.y, -1.0);
	EXPECT_EQ(p1.force.z, -1.0);

	EXPECT_EQ(p2.force.x,  1.0);
	EXPECT_EQ(p2.force.y,  1.0);
	EXPECT_EQ(p2.force.z,  1.0);
}



