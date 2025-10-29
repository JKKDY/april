#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;

using PID = env::internal::ParticleID;

// simple helper to make a dummy particle
static env::internal::Particle make_alive_particle() {
	env::internal::Particle p;
	p.id = 0;
	p.position = {5.0, 5.0, 5.0};
	p.velocity = {0,0,0};
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

// Direct application should mark particle DEAD
TEST(AbsorbBoundaryTest, Apply_SetsParticleDead) {
	const Absorb absorb;
	const env::Box box {{0,0,0}, {10,10,10}};

	auto p = make_alive_particle();
	absorb.apply(p, box, Face::XPlus);

	EXPECT_EQ(p.state, env::ParticleState::DEAD)
		<< "Absorb boundary should mark particle as DEAD";
}

// Topology sanity: outside region, not coupled, no force wrap, no position change
TEST(AbsorbBoundaryTest, Topology_IsOutsideAndPassive) {
	const Absorb absorb;
	const boundary::Topology& topo = absorb.topology;

	EXPECT_LT(topo.boundary_thickness, 0.0)
		<< "Absorb boundaries should have negative thickness (outside domain)";
	EXPECT_FALSE(topo.couples_axis);
	EXPECT_FALSE(topo.force_wrap);
	EXPECT_FALSE(topo.may_change_particle_position);
}


TEST(AbsorbBoundaryTest, CompiledBoundary_Apply_SetsParticleDead) {
	std::variant<Absorb> absorb = Absorb();
	env::Domain domain{{0,0,0}, {10,10,10}};

	// Compile boundary for X+ face
	auto compiled = boundary::internal::compile_boundary(absorb, env::Box::from_domain(domain), Face::XPlus);

	auto p = make_alive_particle();
	env::Box box{{0,0,0}, {10,10,10}};

	compiled.apply(p, box, Face::XPlus);

	EXPECT_EQ(p.state, env::ParticleState::DEAD);
}



template <class ContainerT>
class AbsorbBoundarySystemTestT : public testing::Test {};

using ContainerTypes = testing::Types<DirectSum, LinkedCells>;
TYPED_TEST_SUITE(AbsorbBoundarySystemTestT, ContainerTypes);


// Particle inside the domain should remain alive
TYPED_TEST(AbsorbBoundarySystemTestT, InsideDomain_RemainsAlive) {
	Environment env(forces<NoForce>, boundary::boundaries<Absorb>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Particle in center of domain
	env.add_particle({.id=0, .type=0, .position={5,5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE});

	// Set Absorb on all faces
	env.set_boundaries(Absorb(), all_faces);


	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id0 = mappings.id_map.at(0);
	const auto& p = sys.get_particle_by_index(id0);

	EXPECT_EQ(p.state, ParticleState::ALIVE)
		<< "Particle inside domain should not be affected by absorbing boundaries.";
}

// Absorb boundary full-pipeline test: one particle per face
TYPED_TEST(AbsorbBoundarySystemTestT, EachFace_ParticleMarkedDead) {
	Environment env(forces<NoForce>, boundary::boundaries<Absorb>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Place one particle near each face, moving outward
	env.add_particle({.id=0, .type=0, .position={0.4,5,5},  .velocity={-1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X−
	env.add_particle({.id=1, .type=0, .position={9.6,5,5},  .velocity={+1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X+
	env.add_particle({.id=2, .type=0, .position={5,0.4,5},  .velocity={0,-1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y−
	env.add_particle({.id=3, .type=0, .position={5,9.6,5},  .velocity={0,+1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y+
	env.add_particle({.id=4, .type=0, .position={5,5,0.4},  .velocity={0,0,-1}, .mass=1, .state=ParticleState::ALIVE}); // Z−
	env.add_particle({.id=5, .type=0, .position={5,5,9.6},  .velocity={0,0,+1}, .mass=1, .state=ParticleState::ALIVE}); // Z+

	// Set Absorb on all faces
	env.set_boundaries(Absorb(), all_faces);

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Simulate one time step: move each particle beyond its face
	for (auto pid = sys.index_start(); pid < sys.index_end(); ++pid) {
		auto& p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // crosses boundary
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	// Verify all particles are DEAD
	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.id_map.at(uid);
		const auto& p = sys.get_particle_by_index(iid);
		EXPECT_EQ(p.state, ParticleState::DEAD)
			<< "Particle " << uid << " crossing face " << uid
			<< " should be marked DEAD by Absorb boundary.";
	}
}