#include <gtest/gtest.h>
#include <gmock/gmock.h>


using testing::AnyOf;
using testing::Eq;

#include "april/april.h"
#include "orbit_monitor.h"
#include "constant_force.h"
#include "utils.h"

// TODO tests for get_particle_by_id
using namespace april;


TEST(LinkedCellsTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
	e.add_particle(make_particle(0, {1,2,3}, {}, 1, ParticleState::ALIVE, 0));
	e.add_force(NoForce(), to_type(0));
	e.set_extent({4,4,4});

	auto sys = build_system(e, LinkedCells(4));
    sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_SameCell) {
    Environment e(forces<ConstantForce>);
	e.set_extent({2,2,2});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1.5,0,0}, {}, 2, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(2));
	sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,1});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1.5,0,0}, {}, 2, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(1));
	sys.update_forces();

    auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NoNeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,0.5});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1,0,0}, {}, 1, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(0.5));
	sys.update_forces();
    auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_EQ(p1.force, vec3(0,0,0));
	EXPECT_EQ(p2.force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_IdSpecificForce) {
    Environment e(forces<NoForce, ConstantForce>);
	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 42));
	e.add_particle(make_particle(0, {0,1,0}, {}, 1, ParticleState::ALIVE, 99));
	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));
	e.auto_domain(2);

	auto sys = build_system(e, LinkedCells());
	sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TEST(LinkedCellsTest, TwoParticles_InverseSquare) {
    Environment e(forces<NoForce, Gravity>);

    e.set_extent({10,10,10});

	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(1, {2,0,0}, {}, 2, ParticleState::ALIVE, 1));


	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));

	e.add_force(Gravity(5.0), between_types(0, 1));

	auto sys = build_system(e, LinkedCells());
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


TEST(LinkedCellsTest, OrbitTest) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env (forces<Gravity>);

	env.add_particle(make_particle(0, {0,R,0}, {v, 0, 0}, m));
	env.add_particle(make_particle(0, {0,0,0}, {0, 0, 0}, M));

	env.add_force(Gravity(G), to_type(0));

	env.set_origin({-1.5*v,-1.5*v,0});
	env.set_extent({3*v,3*v,1});

	auto sys = build_system(env, LinkedCells(v));
	sys.update_forces();

	VelocityVerlet integrator(sys, monitor::monitors<OrbitMonitor>);
	integrator.add_monitor(OrbitMonitor(v, R));
	integrator.run_for_duration(0.001, T);

	auto particles = export_particles(sys);

	auto p1 =  particles[0].mass == m ? particles[0] : particles[1];
	auto p2 =  particles[0].mass == M ? particles[0] : particles[1];

	EXPECT_NEAR(p1.velocity.norm(), v, 1e-3);

	EXPECT_NEAR(p1.position.x, 0, 1e-3);
	EXPECT_NEAR(p1.position.y, R, 1e-3);
	EXPECT_EQ(p1.position.z, 0);

	EXPECT_NEAR(p1.velocity.x, v, 1e-3);
	EXPECT_NEAR(p1.velocity.y, 0, 1e-3);
	EXPECT_EQ(p1.velocity.z, 0);

	EXPECT_NEAR(p2.position.x, 0, 1e-3);
	EXPECT_NEAR(p2.position.y, 0, 1e-3);
	EXPECT_NEAR(p2.position.z, 0, 1e-3);

	EXPECT_NEAR(p2.velocity.x, 0, 1e-3);
	EXPECT_NEAR(p2.velocity.y, 0, 1e-3);
	EXPECT_NEAR(p2.velocity.z, 0, 1e-3);
}


TEST(LinkedCellsTest, CollectIndicesInRegion) {
	// Create a simple 3x3x3 grid of particles (27 total)
	auto cuboid = ParticleCuboid{}
		.at(vec3(0.25))
		.velocity({0, 0, 0})
		.count({3, 3, 3})
		.mass(1.0)
		.spacing(1)
		.type(0);

    // Loop through different cell size hints to verify consistency
    for (double cell_size : {0.5, 1.0, 2.0, 5.0}) {

        Environment e(forces<NoForce>);
        e.set_origin({0, 0, 0});
        e.set_extent({5, 5, 5});
		e.add_particles(cuboid);
        e.add_force(NoForce(), to_type(0));

        auto sys = build_system(e, LinkedCells(cell_size));

        // Case 1: small inner region (should include one particle)
        {
            env::Domain region({0.1, 0.1, 0.1}, {0.9, 0.9, 0.9});
            auto indices = sys.query_region(region);
            ASSERT_EQ(indices.size(), 1u);
	        auto p = export_particles(sys)[indices[0]];
            EXPECT_EQ(p.position.x, 0.25);
            EXPECT_EQ(p.position.y, 0.25);
            EXPECT_EQ(p.position.z, 0.25);
        }

        // Case 2: mid region (should include all 27)
        {
            env::Domain region({0, 0, 0}, {5, 5, 5});
            auto indices = sys.query_region(region);
            EXPECT_EQ(indices.size(), 27u);
        }

        // Case 3: partially overlapping region
        {
            env::Domain region({1.5, 1.5, 1.5}, {4.5, 4.5, 4.5});
            std::vector indices = sys.query_region(region);
            EXPECT_GT(indices.size(), 0u);
            EXPECT_LT(indices.size(), 27u);

        	std::unordered_set inside (indices.begin(), indices.end());

			for (size_t id = sys.min_id(); id < sys.max_id(); id++) {
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
            auto indices = sys.query_region(region);
            EXPECT_TRUE(indices.empty());
        }
    }
}

// does nothing except signaling the container to be periodic
struct DummyPeriodicBoundary final : Boundary {
	static constexpr env::FieldMask fields = +env::Field::none;

	DummyPeriodicBoundary()
	: Boundary(0.0, false, true, false ) {}

	template<env::FieldMask M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> &, const env::Box &, const Face) const noexcept{
	}
};

TEST(LinkedCellsTest, PeriodicForceWrap_X) {
	// Iterate over several cell sizes (smaller, medium, larger than extent/2)
	for (double cell_size_hint : {1.0, 3.3, 9.9}) {
		Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
		e.set_origin({0,0,0});
		e.set_extent({10,10,10}); // domain box 10x10x10

		// Two particles, near opposite faces along x
		e.add_particle(make_particle(0, {0.5, 5, 5}, {}, 1, ParticleState::ALIVE, 0));
		e.add_particle(make_particle(0, {9.5, 5, 5}, {}, 1, ParticleState::ALIVE, 1));

		// Simple harmonic force
		e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

		// Enable periodic boundaries on both x faces
		e.set_boundaries(DummyPeriodicBoundary(), {Face::XMinus, Face::XPlus});

		BuildInfo mapping;
		auto sys = build_system(e, LinkedCells(cell_size_hint), &mapping);
		sys.update_forces();

		auto const& out = export_particles(sys);
		ASSERT_EQ(out.size(), 2u);

		auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
		auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

		// They should feel equal and opposite forces due to wrapping
		EXPECT_EQ(p1.force, -p2.force);

		EXPECT_NEAR(p1.force.x, 1.0, 1e-12);
		EXPECT_NEAR(p2.force.x, -1.0, 1e-12);
	}
}


TEST(LinkedCellsTest, PeriodicForceWrap_AllAxes) {
	for (double cell_size_hint : {1.0, 3.3, 9.9}) {

		Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
		e.set_origin({0,0,0});
		e.set_extent({10,10,10});

		// Particles at opposite corners
		e.add_particle(make_particle(0, {0.5, 0.5, 0.5}, {}, 1, ParticleState::ALIVE, 0));
		e.add_particle(make_particle(0, {9.5, 9.5, 9.5}, {}, 1, ParticleState::ALIVE, 1));

		e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

		// Enable full periodicity on all faces
		e.set_boundaries(DummyPeriodicBoundary(), {
			Face::XMinus, Face::XPlus,
			Face::YMinus, Face::YPlus,
			Face::ZMinus, Face::ZPlus
		});

		BuildInfo mapping;
		auto sys = build_system(e, LinkedCells(cell_size_hint), &mapping);
		sys.update_forces();

		auto const& out = export_particles(sys);
		ASSERT_EQ(out.size(), 2u);

		auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
		auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

		// Forces must be equal and opposite
		EXPECT_EQ(p1.force, -p2.force);

		EXPECT_NEAR(p1.force.x, 1.0, 1e-12);
		EXPECT_NEAR(p1.force.y, 1.0, 1e-12);
		EXPECT_NEAR(p1.force.z, 1.0, 1e-12);

		EXPECT_NEAR(p2.force.x, -1.0, 1e-12);
		EXPECT_NEAR(p2.force.y, -1.0, 1e-12);
		EXPECT_NEAR(p2.force.z, -1.0, 1e-12);
	}
}

