#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;

using PID = env::internal::ParticleID;

inline env::internal::Particle make_particle(const vec3& pos, const vec3& vel = {0,0,0}) {
	env::internal::Particle p;
	p.id = 0;
	p.position = pos + vel;
	p.old_position = pos;
	p.velocity = vel;
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

// Direct Application Tests
TEST(PeriodicBoundaryTest, Apply_WrapsAcrossDomain_XPlus) {
	const Periodic periodic;
	const env::Box box({0,0,0}, {10,10,10});

	// Particle just beyond +X boundary
	auto p = make_particle({10.2, 5.0, 5.0});
	periodic.apply(p, box, Face::XPlus);

	EXPECT_NEAR(p.position.x, 0.2, 1e-12);
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}

TEST(PeriodicBoundaryTest, Apply_WrapsAcrossDomain_XMinus) {
	const Periodic periodic;
	const env::Box box({0,0,0}, {10,10,10});

	// Particle just beyond -X boundary
	auto p = make_particle({-0.3, 5.0, 5.0});
	periodic.apply(p, box, Face::XMinus);

	EXPECT_NEAR(p.position.x, 9.7, 1e-12);
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}

TEST(PeriodicBoundaryTest, Apply_WrapsEachAxisCorrectly) {
	const Periodic periodic;
	const env::Box box({0,0,0}, {10,10,10});

	const std::array start_positions = {
		vec3{-0.1, 5, 5},   // X-
		vec3{10.1, 5, 5},   // X+
		vec3{5, -0.2, 5},   // Y-
		vec3{5, 10.3, 5},   // Y+
		vec3{5, 5, -0.4},   // Z-
		vec3{5, 5, 10.5}    // Z+
	};

	constexpr std::array faces = {
		Face::XMinus, Face::XPlus,
		Face::YMinus, Face::YPlus,
		Face::ZMinus, Face::ZPlus
	};

	const std::array expected = {
		vec3{9.9,5,5}, vec3{0.1,5,5},
		vec3{5,9.8,5}, vec3{5,0.3,5},
		vec3{5,5,9.6}, vec3{5,5,0.5}
	};

	for (size_t i = 0; i < faces.size(); ++i) {
		auto p = make_particle(start_positions[i]);
		periodic.apply(p, box, faces[i]);
		EXPECT_NEAR(p.position.x, expected[i].x, 1e-12);
		EXPECT_NEAR(p.position.y, expected[i].y, 1e-12);
		EXPECT_NEAR(p.position.z, expected[i].z, 1e-12);
	}
}

// Topology Sanity
TEST(PeriodicBoundaryTest, Topology_IsOutsideCoupledAndWrapsForces) {
	const Periodic periodic;
	const auto& topo = periodic.topology;

	EXPECT_LT(topo.boundary_thickness, 0.0)
		<< "Periodic boundaries operate outside the domain (teleport wrap).";
	EXPECT_TRUE(topo.couples_axis)
		<< "Periodic boundaries couple opposite faces.";
	EXPECT_TRUE(topo.force_wrap)
		<< "Periodic boundaries enable container force wrapping.";
	EXPECT_TRUE(topo.may_change_particle_position)
		<< "Periodic boundaries may adjust particle positions (teleport).";
}

// 3. Compiled Boundary Variant
TEST(PeriodicBoundaryTest, CompiledBoundary_Apply_WrapsCorrectly) {
	std::variant<Periodic> variant = Periodic();
	env::Domain domain({0,0,0}, {10,10,10});

	auto compiled = boundary::internal::compile_boundary(variant, domain, Face::ZPlus);
	env::internal::Particle p = make_particle({5,5,10.2});
	env::Box box({0,0,0}, {10,10,10});

	compiled.apply(p, box, Face::ZPlus);

	EXPECT_NEAR(p.position.z, 0.2, 1e-12)
		<< "Periodic boundary should wrap Z+ back into domain.";
	EXPECT_NEAR(p.position.x, 5.0, 1e-12);
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
}

// System-level Pipeline Tests
template <class ContainerT>
class PeriodicBoundarySystemTestT : public testing::Test {};

using ContainerTypes = testing::Types<DirectSum, LinkedCells>;
TYPED_TEST_SUITE(PeriodicBoundarySystemTestT, ContainerTypes);

TYPED_TEST(PeriodicBoundarySystemTestT, EachFace_WrapsPositionsAcrossDomain) {
	Environment env(forces<NoForce>, boundary::boundaries<Periodic>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// One particle near each face moving outward
	env.add_particle({.id=0, .type=0, .position={0.4,5,5},  .velocity={-1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X−
	env.add_particle({.id=1, .type=0, .position={9.6,5,5},  .velocity={+1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X+
	env.add_particle({.id=2, .type=0, .position={5,0.4,5},  .velocity={0,-1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y−
	env.add_particle({.id=3, .type=0, .position={5,9.6,5},  .velocity={0,+1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y+
	env.add_particle({.id=4, .type=0, .position={5,5,0.4},  .velocity={0,0,-1}, .mass=1, .state=ParticleState::ALIVE}); // Z−
	env.add_particle({.id=5, .type=0, .position={5,5,9.6},  .velocity={0,0,+1}, .mass=1, .state=ParticleState::ALIVE}); // Z+

	// Enable periodic boundaries on all faces
	env.set_boundaries(Periodic(), all_faces);

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Simulate one integration step: move each particle outside its face
	for (auto pid = sys.index_start(); pid < sys.index_end(); ++pid) {
		auto& p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity;
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	// Expected positions after wrapping (10x10x10 domain)
	std::array expected = {
		vec3{9.4, 5, 5},  // X− → reappears near +X side
		vec3{0.6, 5, 5},  // X+ → reappears near −X side
		vec3{5, 9.4, 5},  // Y− → reappears near +Y side
		vec3{5, 0.6, 5},  // Y+ → reappears near −Y side
		vec3{5, 5, 9.4},  // Z− → reappears near +Z side
		vec3{5, 5, 0.6}   // Z+ → reappears near −Z side
	};

	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.usr_ids_to_impl_ids.at(uid);
		const auto& p = sys.get_particle_by_index(iid);
		EXPECT_NEAR(p.position.x, expected[iid].x, 1e-12);
		EXPECT_NEAR(p.position.y, expected[iid].y, 1e-12);
		EXPECT_NEAR(p.position.z, expected[iid].z, 1e-12);
	}
}


TYPED_TEST(PeriodicBoundarySystemTestT, Integration_CrossAndWrapMaintainsContinuity) {
	Environment env(forces<NoForce>, boundary::boundaries<Periodic>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Single particle heading out +X
	env.add_particle({.id=0, .type=0, .position={9.8,5,5}, .velocity={+1,0,0}, .mass=1, .state=ParticleState::ALIVE});
	env.set_boundaries(Periodic(), all_faces);

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	for (auto pid = sys.index_start(); pid < sys.index_end(); ++pid) {
		auto& p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // simulate step crossing X+
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	const auto& p = sys.get_particle_by_index(mappings.usr_ids_to_impl_ids.at(0));
	EXPECT_NEAR(p.position.x, 0.8, 1e-12)
		<< "Particle crossing +X should reappear at +0.8 within the domain.";
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}
