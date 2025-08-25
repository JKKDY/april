#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <april/env/environment.h>
#include "april/common.h"
#include "april/algo/direct_sum.h"
#include "april/algo/linked_cells.h"
#include "april/core/stoermer_verlet.h"
#include "april/io/output.h"

using namespace april;
using namespace april::env;
using namespace april::core;
using namespace april::io;




TEST(StoermerVerletTest,ConstructionTest) {

	Environment env;
	env.add_particle({}, {}, 1);
	env.add_particle({}, {}, 1);
	env.add_force_to_type(NoForce(), 0);

	constexpr auto algo = algo::DirectSum();
	auto system = compile(env, algo);

	StoermerVerlet<> integrator(system);
	integrator.run_steps(0.1, 10);

	for (auto & p : system.export_particles()) {
		EXPECT_EQ(p.position, vec3(0,0,0));
		EXPECT_EQ(p.velocity, vec3(0,0,0));
	}
}

TEST(StoermerVerletTest, SingleStepNoForceTest) {
	Environment env;
	env.add_particle({}, {1,2,3}, 1);
	env.add_particle({}, {4,5,6}, 2);
	env.add_force_to_type(NoForce(), 0);

	constexpr auto algo = algo::DirectSum();
	auto system = compile(env, algo);

	StoermerVerlet<> integrator(system);
	integrator.run_steps(1, 1);

	std::vector<env::impl::ParticleView> particles;
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


TEST(StoermerVerletTest, SingleStepWithForceTest) {
	Environment env;
	env.add_particle({-1,0,0}, {}, 1 );
	env.add_particle({1,0,0}, {}, 1);
	env.add_force_to_type(InverseSquare(1), 0);

	constexpr auto algo = algo::DirectSum();
	auto system = compile(env, algo);

	StoermerVerlet<> integrator(system);
	integrator.run_steps(0.1, 1);

	std::vector<env::impl::ParticleView> particles;
	for (auto & p : system.export_particles()) {
		particles.push_back(p);
	}

	constexpr double f_mag = 1.0/2/2;

	const auto p1 =  particles[0].position.x < 0 ? particles[0] : particles[1];
	const auto p2 =  particles[0].position.x > 0 ? particles[0] : particles[1];

	EXPECT_EQ(p1.force, vec3(f_mag,0,0));
	EXPECT_EQ(p2.force, vec3(-f_mag,0,0));

	constexpr double vel = 0.1 / 2 * f_mag;

	EXPECT_EQ(p1.velocity, vec3(vel,0,0));
	EXPECT_EQ(p2.velocity, vec3(-vel,0,0));
}


class OrbitMonitor final : public Monitor {
public:
	OrbitMonitor(): Monitor(1) {}
	explicit OrbitMonitor(const double v, const double r): Monitor(1), v(v), r(r) {}

	void record(size_t , double, const Particles& particles) const {
		const auto p = particles[0].mass < 1 ?  particles[0]: particles[1];

		EXPECT_NEAR(p.velocity.norm(), v, 1e-3);
		EXPECT_NEAR(p.position.norm(), r, 1e-3);
	}

	double v{};
	double r{};
};

TEST(StoermerVerletTest, OrbitTest) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env;
	env.add_particle({0,R,0}, {v, 0, 0}, m);
	env.add_particle({0,0,0}, {0, 0, 0}, M);
	env.add_force_to_type(InverseSquare(G), 0);

	constexpr auto algo = algo::DirectSum();
	auto system = compile(env, algo);

	StoermerVerlet<OrbitMonitor> integrator(system);
	integrator.add_monitor(OrbitMonitor(v, R));
	integrator.run(0.001, T);

	std::vector<env::impl::ParticleView> particles;
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
