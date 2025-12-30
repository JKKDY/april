#include <gtest/gtest.h>
#include "april/april.hpp"
#include "utils.h"

using namespace april;



inline env::internal::ParticleRecord<env::NoUserData> make_particle(const vec3& pos, const vec3& vel = {0,0,0}) {
	env::internal::ParticleRecord<env::NoUserData> p;
	p.id = 0;
	p.position = pos + vel;
	p.old_position = pos;
	p.velocity = vel;
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


// Direct Application Tests
TEST(PeriodicBoundaryTest, Apply_WrapsAcrossDomain_XPlus) {
	const Periodic periodic;
	constexpr env::FieldMask Mask = Periodic::fields;

	const env::Box box({0,0,0}, {10,10,10});

	// Particle just beyond +X boundary
	auto p = make_particle({10.2, 5.0, 5.0});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	periodic.apply(ref, box, Face::XPlus);

	EXPECT_NEAR(p.position.x, 0.2, 1e-12);
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}

TEST(PeriodicBoundaryTest, Apply_WrapsAcrossDomain_XMinus) {
	const Periodic periodic;
	constexpr env::FieldMask Mask = Periodic::fields;
	const env::Box box({0,0,0}, {10,10,10});

	// Particle just beyond -X boundary
	auto p = make_particle({-0.3, 5.0, 5.0});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	periodic.apply(ref, box, Face::XMinus);

	EXPECT_NEAR(p.position.x, 9.7, 1e-12);
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}

TEST(PeriodicBoundaryTest, Apply_WrapsEachAxisCorrectly) {
	const Periodic periodic;
	const env::Box box({0,0,0}, {10,10,10});
	constexpr env::FieldMask Mask = Periodic::fields;


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
		auto src = make_source<Mask>(p);
		env::ParticleRef<Mask, env::NoUserData> ref(src);
		periodic.apply(ref, box, faces[i]);
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
	constexpr env::FieldMask Mask = Periodic::fields;
	env::Domain domain({0,0,0}, {10,10,10});

	auto compiled = boundary::internal::compile_boundary(variant, env::Box::from_domain(domain), Face::ZPlus);

	auto p = make_particle({5,5,10.2});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	env::Box box({0,0,0}, {10,10,10});

	compiled.dispatch([&](auto && bc) {
		bc.apply(ref, box, Face::ZPlus);
	});

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
	env.add_particle(make_particle(0, {0.4,5,5}, {-1,0,0}, 1, ParticleState::ALIVE, 0)); // X−
	env.add_particle(make_particle(0, {9.6,5,5}, {+1,0,0}, 1, ParticleState::ALIVE, 1)); // X+
	env.add_particle(make_particle(0, {5,0.4,5}, {0,-1,0}, 1, ParticleState::ALIVE, 2)); // Y−
	env.add_particle(make_particle(0, {5,9.6,5}, {0,+1,0}, 1, ParticleState::ALIVE, 3)); // Y+
	env.add_particle(make_particle(0, {5,5,0.4}, {0,0,-1}, 1, ParticleState::ALIVE, 4)); // Z−
	env.add_particle(make_particle(0, {5,5,9.6}, {0,0,+1}, 1, ParticleState::ALIVE, 5)); // Z+

	// Enable periodic boundaries on all faces
	env.set_boundaries(Periodic(), all_faces);

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Simulate one integration step: move each particle outside its face
	simulate_single_step(sys);

	sys.rebuild_structure();
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
		auto iid = mappings.id_map.at(uid);
		const auto p = get_particle_by_id(sys, iid);

		EXPECT_NEAR(p.position.x, expected[uid].x, 1e-12);
		EXPECT_NEAR(p.position.y, expected[uid].y, 1e-12);
		EXPECT_NEAR(p.position.z, expected[uid].z, 1e-12);
	}
}


TYPED_TEST(PeriodicBoundarySystemTestT, Integration_CrossAndWrapMaintainsContinuity) {
	Environment env(forces<NoForce>, boundary::boundaries<Periodic>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Single particle heading out +X
	env.add_particle(make_particle(0, {9.8,5,5}, {+1,0,0}, 1, ParticleState::ALIVE, 0));
	env.set_boundaries(Periodic(), all_faces);

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	simulate_single_step(sys);


	sys.rebuild_structure();
	sys.apply_boundary_conditions();

	const auto p = get_particle(sys, mappings.id_map.at(0));

	EXPECT_NEAR(p.position.x, 0.8, 1e-12)
		<< "Particle crossing +X should reappear at +0.8 within the domain.";
	EXPECT_NEAR(p.position.y, 5.0, 1e-12);
	EXPECT_NEAR(p.position.z, 5.0, 1e-12);
}
