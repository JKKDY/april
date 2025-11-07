#include <gtest/gtest.h>
#include "april/april.h"
#include "utils.h"

using namespace april;

// simple helper to make a dummy particle
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

// Direct application should reflect the particles position
TEST(ReflectiveBoundaryTest, Apply_InvertsVelocityAndReflectsPosition) {
	const Reflective reflective;
	const env::Box box({0,0,0}, {10,10,10});

	// heading out X+. Intersection at {10, 5, 5}
	auto p = make_particle({9.5,4.5,4.5}, {2,2,2});
	auto pf = env::internal::ParticleRecordFetcher(p);
	reflective.apply(pf, box, Face::XPlus);

	EXPECT_TRUE(box.contains(p.position));
	EXPECT_EQ(p.position.x, 8.5);
	EXPECT_EQ(p.position.y, 6.5);
	EXPECT_EQ(p.position.z, 6.5);
	EXPECT_EQ(p.velocity.x, -2.0);
	EXPECT_EQ(p.velocity.y, 2.0);
	EXPECT_EQ(p.velocity.z, 2.0);
}


// Topology sanity: outside region, not coupled, no force wrap, position change
TEST(ReflectiveBoundaryTest, Topology_IsOutsideAndChangesPosition) {
	const Reflective reflective;
	const auto& topology = reflective.topology;

	EXPECT_LT(topology.boundary_thickness, 0.0)
		<< "Reflective boundaries operate outside the domain (negative thickness).";
	EXPECT_FALSE(topology.couples_axis);
	EXPECT_FALSE(topology.force_wrap);
	EXPECT_TRUE(topology.may_change_particle_position)
		<< "Reflective boundary should adjust particle position.";
}


TEST(AbsorbBoundaryTest, CompiledBoundary_Apply_InvertsVelocityAndReflectsPosition) {
	std::variant<Reflective> reflect = Reflective();
	env::Domain domain({0,0,0}, {10,10,10});

	// Compile boundary for X+ face
	auto compiled = boundary::internal::compile_boundary(reflect, env::Box::from_domain(domain), Face::XPlus);

	auto p = make_particle({9.8,5,5}, {+1,0,0});
	env::Box box{{0,0,0}, {10,10,10}};

	auto pf = env::internal::ParticleRecordFetcher(p);
	compiled.apply(pf, box, Face::XPlus);

	EXPECT_TRUE(box.contains(p.position));
	EXPECT_EQ(p.position.x, 9.2);
	EXPECT_EQ(p.position.y, 5);
	EXPECT_EQ(p.position.z, 5);
	EXPECT_NEAR(p.velocity.x, -1.0, 1e-12);
	EXPECT_NEAR(p.velocity.y, 0.0, 1e-12);
	EXPECT_NEAR(p.velocity.z, 0.0, 1e-12);
}



template <class ContainerT>
class ReflectiveBoundarySystemTestT : public testing::Test {};

using ContainerTypes = testing::Types<DirectSum, LinkedCells>;
TYPED_TEST_SUITE(ReflectiveBoundarySystemTestT, ContainerTypes);


TYPED_TEST(ReflectiveBoundarySystemTestT, EachFace_ReflectsVelocityInNormal) {
    Environment env(forces<NoForce>, boundary::boundaries<Reflective>);
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

	env.set_boundaries(Reflective(), all_faces);


    BuildInfo mappings;
    auto sys = build_system(env, TypeParam(), &mappings);

    // simulate one step
	simulate_single_step(sys);


    sys.register_all_particle_movements();
    sys.apply_boundary_conditions();

	// expected positions
	std::array expected_pos = {
		vec3{0.6, 5, 5},
		vec3{9.4, 5, 5},
		vec3{5, 0.6, 5},
		vec3{5, 9.4, 5},
		vec3{5, 5, 0.6},
		vec3{5, 5, 9.4},
	};

    // Expected velocity signs after reflection
    std::array expected_vel = {
        vec3{+1,0,0},
    	vec3{-1,0,0},
        vec3{0,+1,0},
    	vec3{0,-1,0},
        vec3{0,0,+1},
    	vec3{0,0,-1}
    };

    for (int uid = 0; uid < 6; ++uid) {
        auto iid = mappings.id_map.at(uid);

    	const auto p = get_particle(sys, iid);

    	EXPECT_EQ(p.position.x, expected_pos[iid].x);
    	EXPECT_EQ(p.position.y, expected_pos[iid].y);
    	EXPECT_EQ(p.position.z, expected_pos[iid].z);

        EXPECT_EQ(p.velocity.x, expected_vel[iid].x);
        EXPECT_EQ(p.velocity.y, expected_vel[iid].y);
        EXPECT_EQ(p.velocity.z, expected_vel[iid].z);
    }
}