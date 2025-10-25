#include <gtest/gtest.h>
#include <gmock/gmock.h>


#include "april/april.h"
#include "OrbitMonitor.h"

using namespace april;

TEST(Yoshida4Test,ConstructionTest) {

	Environment env (forces<NoForce>);
	env.add_particle({}, {}, 1);
	env.add_particle({}, {}, 1);
	env.add_force(NoForce(), to_type(0));
	env.set_extent({2,2,2});
	env.set_origin({-1,-1,-1});

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	Yoshida4 integrator(system);
	integrator.run_steps(0.1, 10);

	for (auto & p : system.export_particles()) {
		EXPECT_EQ(p.position, vec3(0,0,0));
		EXPECT_EQ(p.velocity, vec3(0,0,0));
	}
}

TEST(Yoshida4Test, SingleStepNoForceTest) {
	Environment env (forces<NoForce>);
	env.add_particle({}, {1,2,3}, 1);
	env.add_particle({}, {4,5,6}, 2);
	env.add_force(NoForce(), to_type(0));
	env.set_extent(20 * vec3(1));
	env.set_origin(10*vec3{-1});

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	Yoshida4 integrator(system);
	integrator.run_steps(1, 1);

	std::vector<ParticleView> particles;
	for (auto & p : system.export_particles()) {
		particles.push_back(p);
	}

	auto p1 =  particles[0].mass == 1 ? particles[0] : particles[1];
	auto p2 =  particles[0].mass == 2 ? particles[0] : particles[1];

	EXPECT_NEAR(p1.position.x, 1, 1e-5);
	EXPECT_NEAR(p1.position.y, 2, 1e-5);
	EXPECT_NEAR(p1.position.z, 3, 1e-5);

	EXPECT_NEAR(p2.position.x, 4, 1e-5);
	EXPECT_NEAR(p2.position.y, 5, 1e-5);
	EXPECT_NEAR(p2.position.z, 6, 1e-5);

	EXPECT_EQ(p1.velocity, vec3(1,2,3));
	EXPECT_EQ(p2.velocity, vec3(4,5,6));
}


TEST(Yoshida4Test, SingleStepWithForceTest) {
	Environment env (forces<PowerLaw>);
	env.add_particle({-1,0,0}, {}, 1 );
	env.add_particle({1,0,0}, {}, 1);
	env.add_force(PowerLaw(2, 1), to_type(0));
	env.set_extent({4,4,4});
	env.set_origin({-2,-2,-2});

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	Yoshida4 integrator(system);
	integrator.run_steps(0.1, 1);

	std::vector<ParticleView> particles;
	for (auto & p : system.export_particles()) {
		particles.push_back(p);
	}

	constexpr double f_mag = 1.0/2/2;

	const auto p1 =  particles[0].position.x < 0 ? particles[0] : particles[1];
	const auto p2 =  particles[0].position.x > 0 ? particles[0] : particles[1];

	EXPECT_NEAR(p1.force.x,  f_mag, 1e-2);
	EXPECT_NEAR(p1.force.y,  0.0,   1e-2);
	EXPECT_NEAR(p1.force.z,  0.0,   1e-2);

	EXPECT_NEAR(p2.force.x, -f_mag, 1e-2);
	EXPECT_NEAR(p2.force.y,  0.0,   1e-2);
	EXPECT_NEAR(p2.force.z,  0.0,   1e-2);

	constexpr double vel = 0.1 / 2 * f_mag;

	EXPECT_NEAR(p1.velocity.x,  vel, 1e-2);
	EXPECT_NEAR(p1.velocity.y,  0.0, 1e-2);
	EXPECT_NEAR(p1.velocity.z,  0.0, 1e-2);

	EXPECT_NEAR(p2.velocity.x, -vel, 1e-2);
	EXPECT_NEAR(p2.velocity.y,  0.0, 1e-2);
	EXPECT_NEAR(p2.velocity.z,  0.0, 1e-2);
}


TEST(Yoshida4Test, OrbitTest) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env (forces<PowerLaw>);
	env.add_particle({0,0,0}, {0, 0, 0}, M);
	env.add_particle({0,R,0}, {v, 0, 0}, m);
	env.add_force(PowerLaw(2, G), to_type(0));
	env.set_extent(vec3{R,R,R}*4);
	env.set_origin(vec3{-R,-R,-R} * 2);

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	Yoshida4 integrator(system, monitors<OrbitMonitor>);
	integrator.add_monitor(OrbitMonitor(v, R));
	integrator.run_for(0.001, T);

	std::vector<ParticleView> particles;
	for (auto & p : system.export_particles()) {
		particles.push_back(p);
	}

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


TEST(Yoshida4Test, OrbitTestSplitRuns) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env (forces<PowerLaw>);
	env.add_particle({0,0,0}, {0, 0, 0}, M);
	env.add_particle({0,R,0}, {v, 0, 0}, m);
	env.add_force(PowerLaw(2, G), to_type(0));
	env.set_extent(vec3{R,R,R}*4);
	env.set_origin(vec3{-R,-R,-R} * 2);

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	{
		Yoshida4 integrator(system, monitor::monitors<OrbitMonitor>);
		integrator.add_monitor(OrbitMonitor(v, R));
		integrator.run_for(0.001, T/2);
		EXPECT_NEAR(system.time(), T/2, 0.005);
	}

	{
		Yoshida4 integrator(system, monitor::monitors<OrbitMonitor>);
		integrator.add_monitor(OrbitMonitor(v, R));
		integrator.run_for(0.001, T/2);
		EXPECT_NEAR(system.time(), T, 0.005);
	}

	std::vector<ParticleView> particles;
	for (auto & p : system.export_particles()) {
		particles.push_back(p);
	}

	auto p1 =  particles[0].mass == m ? particles[0] : particles[1];
	auto p2 =  particles[0].mass == M ? particles[0] : particles[1];

	EXPECT_NEAR(p1.velocity.norm(), v, 1e-3);

	// more relaxed conditions since #integration steps may be off by a little
	EXPECT_NEAR(p1.position.x, 0, 2e-3);
	EXPECT_NEAR(p1.position.y, R, 2e-3);
	EXPECT_EQ(p1.position.z, 0);

	EXPECT_NEAR(p1.velocity.x, v, 2e-3);
	EXPECT_NEAR(p1.velocity.y, 0, 2e-3);
	EXPECT_EQ(p1.velocity.z, 0);

	EXPECT_NEAR(p2.position.x, 0, 2e-3);
	EXPECT_NEAR(p2.position.y, 0, 2e-3);
	EXPECT_NEAR(p2.position.z, 0, 2e-3);

	EXPECT_NEAR(p2.velocity.x, 0, 2e-3);
	EXPECT_NEAR(p2.velocity.y, 0, 2e-3);
	EXPECT_NEAR(p2.velocity.z, 0, 2e-3);
}