#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;

using PID = env::internal::ParticleID;

// Dummy test force. always returns constant value within cutoff
struct ConstantForce {
	double value;
	double rc;
	explicit ConstantForce(const double value = 1.0, const double rc = 10.0)
		: value(value), rc(rc) {}

	[[nodiscard]] double cutoff() const noexcept { return rc; }
	[[nodiscard]] double apply(const env::internal::Particle&, const double dist) const noexcept {
		return (dist <= rc) ? value : 0.0;
	}
};


inline env::internal::Particle make_particle(const vec3& pos) {
	env::internal::Particle p;
	p.id = 0;
	p.position = pos;
	p.velocity = {0,0,0};
	p.force = {0,0,0};
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

// Direct application test should add constant force
TEST(RepulsiveBoundaryTest, Apply_AddsInwardForce) {
	ConstantForce f{5.0, 10.0};
	const Repulsive rep(f);

	env::internal::Particle p = make_particle({9.5,5,5});
	const env::Box box({0,0,0}, {10,10,10});

	rep.apply(p, box, Face::XPlus);

	EXPECT_NEAR(p.force.x, -5.0, 1e-12);
	EXPECT_NEAR(p.force.y,  0.0, 1e-12);
	EXPECT_NEAR(p.force.z,  0.0, 1e-12);

	// X- face should push opposite direction
	p.force = {0,0,0};
	rep.apply(p, box, Face::XMinus);
	EXPECT_NEAR(p.force.x, +5.0, 1e-12);
}

// Topology sanity: inside region, not coupled, no force wrap, position change
TEST(RepulsiveBoundaryTest, Topology_IsInsideAndChangesPosition) {
	ConstantForce f{1.0, 3.0};
	const Repulsive rep(f);

	const auto& topology = rep.topology;
	EXPECT_GT(topology.boundary_thickness, 0.0)
		<< "Repulsive boundaries operate inside the domain (positive thickness).";
	EXPECT_FALSE(topology.couples_axis);
	EXPECT_FALSE(topology.force_wrap);
	EXPECT_TRUE(topology.may_change_particle_position);
}

// CompiledBoundary apply test
TEST(RepulsiveBoundaryTest, CompiledBoundary_Apply_AddsInwardForce) {
	ConstantForce f{2.0, 5.0};
	std::variant<std::monostate, Repulsive<ConstantForce>> variant = Repulsive(f);
	env::Domain domain({0,0,0}, {10,10,10});

	auto compiled = boundary::internal::compile_boundary(variant, domain, Face::YMinus);

	env::internal::Particle p = make_particle({5,0.3,5});
	const env::Box box({0,0,0}, {10,10,10});

	compiled.apply(p, box, Face::YMinus);

	EXPECT_NEAR(p.force.y, +2.0, 1e-12)
		<< "Force on Y- face should push inward (+Y direction).";
	EXPECT_NEAR(p.force.x, 0.0, 1e-12);
	EXPECT_NEAR(p.force.z, 0.0, 1e-12);
}


// System-level pipeline test
template <class ContainerT>
class RepulsiveBoundarySystemTestT : public testing::Test {};
using ContainerTypes = testing::Types<DirectSum, LinkedCells>;
TYPED_TEST_SUITE(RepulsiveBoundarySystemTestT, ContainerTypes);


TYPED_TEST(RepulsiveBoundarySystemTestT, EachFace_AppliesInwardForce) {
	ConstantForce f{3.0, 5.0};
	Environment env(forces<NoForce>, boundary::boundaries<Repulsive<ConstantForce>>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Particles near each face
	env.add({.id=0, .type=0, .position={0.5,5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // X-
	env.add({.id=1, .type=0, .position={9.5,5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // X+
	env.add({.id=2, .type=0, .position={5,0.5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Y-
	env.add({.id=3, .type=0, .position={5,9.5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Y+
	env.add({.id=4, .type=0, .position={5,5,0.5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Z-
	env.add({.id=5, .type=0, .position={5,5,9.5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Z+

	env.set_boundaries(std::array{
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f)
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Apply boundaries
	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	// Expected inward forces per face
	const std::array expected = {
		vec3{+3,0,0}, vec3{-3,0,0},
		vec3{0,+3,0}, vec3{0,-3,0},
		vec3{0,0,+3}, vec3{0,0,-3}
	};

	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.usr_ids_to_impl_ids.at(uid);
		const auto& p = sys.get_particle_by_index(iid);
		EXPECT_NEAR(p.force.x, expected[iid].x, 1e-12);
		EXPECT_NEAR(p.force.y, expected[iid].y, 1e-12);
		EXPECT_NEAR(p.force.z, expected[iid].z, 1e-12);
	}
}