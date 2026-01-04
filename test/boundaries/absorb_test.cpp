#include <gtest/gtest.h>
#include "april/april.hpp"

#include "utils.h"
using namespace april;



using PID = ParticleID;

// simple helper to make a dummy particle
static env::internal::ParticleRecord<env::NoUserData> make_alive_particle() {
	env::internal::ParticleRecord<env::NoUserData> p;
	p.id = 0;
	p.position = {5.0, 5.0, 5.0};
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

// Direct application should mark particle DEAD
TEST(AbsorbBoundaryTest, Apply_SetsParticleDead) {
	const Absorb absorb;
	constexpr env::FieldMask Mask = Absorb::fields;

	const env::Box box {{0,0,0}, {10,10,10}};

	auto p = make_alive_particle();
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	absorb.apply(ref, box, Face::XPlus);

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
	constexpr env::FieldMask Mask = Absorb::fields;

	env::Domain domain{{0,0,0}, {10,10,10}};

	// Compile boundary for X+ face
	auto compiled = boundary::internal::compile_boundary(absorb, env::Box::from_domain(domain), Face::XPlus);

	auto p = make_alive_particle();
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	env::Box box{{0,0,0}, {10,10,10}};

	compiled.dispatch([&](auto && bc) {
		bc.apply(ref, box, Face::XPlus);
	});

	EXPECT_EQ(p.state, env::ParticleState::DEAD);
}



template <class ContainerT>
class AbsorbBoundarySystemTestT : public testing::Test {};

using ContainerTypes = testing::Types<DirectSumAoS, LinkedCellsAoS>;
TYPED_TEST_SUITE(AbsorbBoundarySystemTestT, ContainerTypes);


// Particle inside the domain should remain alive
TYPED_TEST(AbsorbBoundarySystemTestT, InsideDomain_RemainsAlive) {
	Environment env(forces<NoForce>, boundary::boundaries<Absorb>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Particle in center of domain
	env.add_particle(make_particle(0, {5,5,5}, {}, 1, ParticleState::ALIVE, 0));

	// Set Absorb on all faces
	env.set_boundaries(Absorb(), all_faces);


	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	sys.rebuild_structure();
	sys.apply_boundary_conditions();

	auto id0 = mappings.id_map.at(0);
	const auto p = get_particle(sys, id0);

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
	env.add_particle(make_particle(0, {0.4,5,5}, {-1,0,0}, 1, ParticleState::ALIVE, 0)); // X−
	env.add_particle(make_particle(0, {9.6,5,5}, {+1,0,0}, 1, ParticleState::ALIVE, 1)); // X+
	env.add_particle(make_particle(0, {5,0.4,5}, {0,-1,0}, 1, ParticleState::ALIVE, 2)); // Y−
	env.add_particle(make_particle(0, {5,9.6,5}, {0,+1,0}, 1, ParticleState::ALIVE, 3)); // Y+
	env.add_particle(make_particle(0, {5,5,0.4}, {0,0,-1}, 1, ParticleState::ALIVE, 4)); // Z−
	env.add_particle(make_particle(0, {5,5,9.6}, {0,0,+1}, 1, ParticleState::ALIVE, 5)); // Z+

	// Set Absorb on all faces
	env.set_boundaries(Absorb(), all_faces);

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Simulate one time step: move each particle beyond its face
	simulate_single_step(sys);

	sys.rebuild_structure();
	sys.apply_boundary_conditions();

	// Verify all particles are DEAD
	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.id_map.at(uid);
		const auto p = get_particle(sys, iid);
		EXPECT_EQ(p.state, ParticleState::DEAD)
			<< "Particle " << uid << " crossing face " << uid
			<< " should be marked DEAD by Absorb boundary.";
	}
}