#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::AnyOf;
using testing::Eq;

#include "../../include/april/april.h"
using namespace april;

// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final {
	vec3 v;
	double cutoff_radius;
	ConstantForce(double x, double y, double z, double cutoff = -1) : v{x,y,z} {
		cutoff_radius = cutoff;
	}
	vec3 operator()(const env::internal::Particle&, const env::internal::Particle&, const vec3&) const noexcept {
		return v;
	}

	[[nodiscard]] ConstantForce mix(const ConstantForce& other) const noexcept {
		return {
			v.x + other.v.x,
			v.y + other.v.y,
			v.z + other.v.z,
			std::max(cutoff_radius, other.cutoff_radius)
		};
	}
};


TEST(LinkedCellsTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
    e.add(Particle{.id = 0, .type = 0, .position={1,2,3},.velocity={0,0,0}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));

	e.set_extent({4,4,4});

	auto sys = build_system(e, LinkedCells(4));
    sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_SameCell) {
    Environment e(forces<ConstantForce>);
	e.set_extent({2,2,2});
	e.set_origin({0,0,0});
	e.add(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 1, .type = 7, .position={1,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(2));
	sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,1});
	e.set_origin({0,0,0});
	e.add(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add(Particle{.id = 1, .type = 7, .position={1.5,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(1));
	sys.update_forces();

	auto const& out = sys.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NoNeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,0.5});
	e.set_origin({0,0,0});
	e.add(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add(Particle{.id = 1, .type = 7, .position={1.5,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, LinkedCells(0.5));
	sys.update_forces();
	auto const& out = sys.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_EQ(p1.force, vec3(0,0,0));
	EXPECT_EQ(p2.force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_IdSpecificForce) {
    Environment e(forces<NoForce, ConstantForce>);
    e.add(Particle{.id = 42, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 99, .type = 0, .position={0,1,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});

	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));

	auto sys = build_system(e, LinkedCells());
	sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TEST(LinkedCellsTest, TwoParticles_InverseSquare) {
    Environment e(forces<NoForce, InverseSquare>);

    e.set_extent({10,10,10});

	e.add(Particle{.id = 0, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 1, .type = 1, .position={2,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});

	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));

	e.add_force(InverseSquare(5.0), between_types(0, 1));

	auto sys = build_system(e, LinkedCells());
	sys.update_forces();

	auto const& out = sys.export_particles();
    // find each
    const auto& pa = (out[0].mass == 1 ? out[0] : out[1]);
    const auto& pb = (out[1].mass == 2 ? out[1] : out[0]);
    // magnitude = pre * m1*m2 / r^3 = 5*1*2/(2^3)=10/8=1.25  direction from pa->pb = (2,0,0)
    // force on pa = 1.25*(2,0,0) = (2.5,0,0); on pb = (-2.5,0,0)
    EXPECT_NEAR(pa.force.x, 2.5, 1e-12);
    EXPECT_NEAR(pb.force.x, -2.5, 1e-12);
    EXPECT_EQ(pa.force.y, 0.0);
    EXPECT_EQ(pb.force.y, 0.0);
}


class OrbitMonitor final : public monitor::Monitor {
public:
	OrbitMonitor(): Monitor(1) {}
	explicit OrbitMonitor(const double v, const double r): Monitor(1), v(v), r(r) {}

	void record(const size_t i, double, const Particles& particles) const {
		const auto p = particles[0].mass < 1 ?  particles[0]: particles[1];

		EXPECT_NEAR(p.velocity.norm(), v, 1e-3) << "Velocity mismatch at step " << i;
		EXPECT_NEAR(p.position.norm(), r, 1e-3) << "Position mismatch at step " << i;

		std::cerr << "step " << i << std::endl;
	}

	double v{};
	double r{};
};

TEST(LinkedCellsTest, OrbitTest) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env (forces<InverseSquare>);
	env.add({0,R,0}, {v, 0, 0}, m);
	env.add({0,0,0}, {0, 0, 0}, M);
	env.add_force(InverseSquare(G), to_type(0));

	env.set_origin({-1.5*v,-1.5*v,0});
	env.set_extent({3*v,3*v,1});

	auto sys = build_system(env, LinkedCells(v));
	sys.update_forces();

	StoermerVerlet integrator(sys, monitor::monitors<OrbitMonitor>);
	integrator.add_monitor(OrbitMonitor(v, R));
	integrator.run_for(0.001, T);

	std::vector<ParticleView> particles;
	for (auto & p : sys.export_particles()) {
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


