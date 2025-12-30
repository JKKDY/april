#include <gtest/gtest.h>
#include "april/april.hpp"

#include "utils.h"

using namespace april;

struct ConstantForce {
	double value;
	double rc;
	explicit ConstantForce(const double value = 1.0, const double rc = 10.0)
		: value(value), rc(rc) {}

	[[nodiscard]] double cutoff() const noexcept { return rc; }

	[[nodiscard]] double apply(const double dist) const noexcept {
		return (dist <= rc) ? value : 0.0;
	}
};

inline env::internal::ParticleRecord<env::NoUserData> make_particle(const vec3& pos) {
	env::internal::ParticleRecord<env::NoUserData> p;
	p.id = 0;
	p.position = pos;
	p.force = {0,0,0};
	p.velocity = {0,0,0};
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

template<env::FieldMask Mask, typename RecordT>
auto make_source(RecordT& record) {
	// Determine constness based on RecordT (allows making const sources from const records)
	constexpr bool IsConst = std::is_const_v<RecordT>;
	using UserDataT = typename RecordT::user_data_t;

	env::ParticleSource<Mask, UserDataT, IsConst> src;

	if constexpr (env::has_field_v<Mask, env::Field::position>)     src.position     = &record.position;
	if constexpr (env::has_field_v<Mask, env::Field::velocity>)     src.velocity     = &record.velocity;
	if constexpr (env::has_field_v<Mask, env::Field::force>)        src.force        = &record.force;
	if constexpr (env::has_field_v<Mask, env::Field::old_position>) src.old_position = &record.old_position;
	if constexpr (env::has_field_v<Mask, env::Field::mass>)         src.mass         = &record.mass;
	if constexpr (env::has_field_v<Mask, env::Field::state>)        src.state        = &record.state;
	if constexpr (env::has_field_v<Mask, env::Field::type>)         src.type         = &record.type;
	if constexpr (env::has_field_v<Mask, env::Field::id>)           src.id           = &record.id;
	if constexpr (env::has_field_v<Mask, env::Field::user_data>)    src.user_data    = &record.user_data;

	return src;
}

// Direct application test should add constant force
TEST(RepulsiveBoundaryTest, Apply_AddsInwardForce) {
	ConstantForce f{5.0, 10.0};
	const Repulsive rep(f);
	constexpr env::FieldMask Mask = Repulsive<ConstantForce>::fields;

	auto p = make_particle({9.5,5,5});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	const env::Box box({0,0,0}, {10,10,10});

	rep.apply(ref, box, Face::XPlus);

	EXPECT_NEAR(p.force.x, -5.0, 1e-12);
	EXPECT_NEAR(p.force.y,  0.0, 1e-12);
	EXPECT_NEAR(p.force.z,  0.0, 1e-12);

	// X- face should push opposite direction
	p.force = {0,0,0};
	rep.apply(ref, box, Face::XMinus);
	EXPECT_NEAR(p.force.x, +5.0, 1e-12);
}

// Topology sanity: inside region, not coupled, no force wrap, np position change
TEST(RepulsiveBoundaryTest, Topology_IsInsideAndNoChangesPosition) {
	ConstantForce f{1.0, 3.0};
	const Repulsive rep(f);

	const auto& topology = rep.topology;
	EXPECT_GT(topology.boundary_thickness, 0.0)
		<< "Repulsive boundaries operate inside the domain (positive thickness).";
	EXPECT_FALSE(topology.couples_axis);
	EXPECT_FALSE(topology.force_wrap);
	EXPECT_FALSE(topology.may_change_particle_position);
}

// CompiledBoundary apply test
TEST(RepulsiveBoundaryTest, CompiledBoundary_Apply_AddsInwardForce) {
	ConstantForce f{2.0, 5.0};
	std::variant<Repulsive<ConstantForce>> variant = Repulsive(f);
	env::Domain domain({0,0,0}, {10,10,10});

	constexpr env::FieldMask Mask = Repulsive<ConstantForce>::fields;

	auto compiled = boundary::internal::compile_boundary(variant, env::Box::from_domain(domain), Face::YMinus);

	auto p = make_particle({5,0.3,5});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	const env::Box box({0,0,0}, {10,10,10});

	compiled.dispatch([&](auto && bc) {
		bc.apply(ref, box, Face::YMinus);
	});

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
	env.add_particle(make_particle(0, {0.5,5,5}, {}, 1, ParticleState::ALIVE, 0)); // X-
	env.add_particle(make_particle(0, {9.5,5,5}, {}, 1, ParticleState::ALIVE, 1)); // X+
	env.add_particle(make_particle(0, {5,0.5,5}, {}, 1, ParticleState::ALIVE, 2)); // Y-
	env.add_particle(make_particle(0, {5,9.5,5}, {}, 1, ParticleState::ALIVE, 3)); // Y+
	env.add_particle(make_particle(0, {5,5,0.5}, {}, 1, ParticleState::ALIVE, 4)); // Z-
	env.add_particle(make_particle(0, {5,5,9.5}, {}, 1, ParticleState::ALIVE, 5)); // Z+

	env.set_boundaries(std::array{
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f)
	});

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Apply boundaries
	sys.rebuild_structure();
	sys.apply_boundary_conditions();

	// Expected inward forces per face
	const std::array expected = {
		vec3{+3,0,0}, vec3{-3,0,0},
		vec3{0,+3,0}, vec3{0,-3,0},
		vec3{0,0,+3}, vec3{0,0,-3}
	};

	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.id_map.at(uid);
		const auto p = get_particle_by_id(sys, iid);
		EXPECT_NEAR(p.force.x, expected[uid].x, 1e-12);
		EXPECT_NEAR(p.force.y, expected[uid].y, 1e-12);
		EXPECT_NEAR(p.force.z, expected[uid].z, 1e-12);
	}
}