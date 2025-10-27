#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;

using PID = env::internal::ParticleID;


// simple helper to make a dummy particle
inline env::internal::Particle make_particle(const vec3& pos, const vec3& vel) {
	env::internal::Particle p;
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
	reflective.apply(p, box, Face::XPlus);

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
	auto compiled = boundary::internal::compile_boundary(reflect, domain, Face::XPlus);

	auto p = make_particle({9.8,5,5}, {+1,0,0});
	env::Box box{{0,0,0}, {10,10,10}};

	compiled.apply(p, box, Face::XPlus);

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
    env.add_particle({.id=0, .type=0, .position={0.4,5,5},  .velocity={-1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X−
    env.add_particle({.id=1, .type=0, .position={9.6,5,5},  .velocity={+1,0,0}, .mass=1, .state=ParticleState::ALIVE}); // X+
    env.add_particle({.id=2, .type=0, .position={5,0.4,5},  .velocity={0,-1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y−
    env.add_particle({.id=3, .type=0, .position={5,9.6,5},  .velocity={0,+1,0}, .mass=1, .state=ParticleState::ALIVE}); // Y+
    env.add_particle({.id=4, .type=0, .position={5,5,0.4},  .velocity={0,0,-1}, .mass=1, .state=ParticleState::ALIVE}); // Z−
    env.add_particle({.id=5, .type=0, .position={5,5,9.6},  .velocity={0,0,+1}, .mass=1, .state=ParticleState::ALIVE}); // Z+

	env.set_boundaries(Reflective(), all_faces);


    UserToInternalMappings mappings;
    auto sys = build_system(env, TypeParam(), &mappings);

    // simulate one step
    for (auto pid = sys.index_start(); pid < sys.index_end(); ++pid) {
        auto& p = sys.get_particle_by_index(pid);
        p.old_position = p.position;
        p.position = p.old_position + p.velocity;
    }

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
        auto iid = mappings.usr_ids_to_impl_ids.at(uid);
    	env::internal::Particle & p = sys.get_particle_by_index(iid);

    	EXPECT_EQ(p.position.x, expected_pos[iid].x);
    	EXPECT_EQ(p.position.y, expected_pos[iid].y);
    	EXPECT_EQ(p.position.z, expected_pos[iid].z);

        EXPECT_EQ(p.velocity.x, expected_vel[iid].x);
        EXPECT_EQ(p.velocity.y, expected_vel[iid].y);
        EXPECT_EQ(p.velocity.z, expected_vel[iid].z);
    }
}