#include <gtest/gtest.h>
#include <april/env/environment.h>
#include "april/common.h"
#include <gmock/gmock.h>

#include "april/containers/linked_cells.h"
#include "april/core/stoermer_verlet.h"
#include "april/io/monitor.h"


using testing::AnyOf;
using testing::Eq;

using namespace april;
using namespace april::env;
using namespace april::core;

struct ConstantForce final : Force {
	vec3 v;
	ConstantForce(double x, double y, double z, double cutoff = -1) : v{x,y,z} {
		cutoff_radius = cutoff;
	}
	vec3 operator()(const env::impl::Particle&, const env::impl::Particle&, const vec3&) const noexcept override {
		return v;
	}
	std::unique_ptr<Force> mix(const Force* other) const override {
		auto* o = dynamic_cast<const ConstantForce*>(other);
		if (!o) throw std::invalid_argument("mix mismatch");
		return std::make_unique<ConstantForce>(
			v.x + o->v.x,
			v.y + o->v.y,
			v.z + o->v.z,
			std::max(cutoff_radius, o->cutoff_radius)
		);
	}
};


TEST(LinkedCellsTest, SingleParticle_NoForce) {
    Environment e;
    e.add_particle(Particle{.id = 0, .type = 0, .position={1,2,3},.velocity={0,0,0}, .mass=1.0, .state=ParticleState::ALIVE});
    e.add_force_to_type(NoForce(), 0);
	e.set_extent({4,4,4});
	e.set_container(std::make_unique<LinkedCells>(4));

    e.build();

    e.update_forces();
    auto const& out = e.export_particles();
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_SameCell) {
    Environment e;
	e.set_extent({2,2,2});
	e.set_origin({0,0,0});
	e.set_container(std::make_unique<LinkedCells>(2));
	e.add_particle(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 1, .type = 7, .position={1,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force_to_type(ConstantForce(3,4,5), 7);
    e.build();

    e.update_forces();
    auto const& out = e.export_particles();
    ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NeighbouringCell) {
	Environment e;
	e.set_extent({2,1,1});
	e.set_origin({0,0,0});
	e.set_container(std::make_unique<LinkedCells>(1));
	e.add_particle(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_particle(Particle{.id = 1, .type = 7, .position={1.5,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force_to_type(ConstantForce(3,4,5), 7);
	e.build();

	e.update_forces();
	auto const& out = e.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, -p2.force);
}

TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NoNeighbouringCell) {
	Environment e;
	e.set_extent({2,1,1});
	e.set_origin({0,0,0});
	e.set_container(std::make_unique<LinkedCells>(0.5));
	e.add_particle(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_particle(Particle{.id = 1, .type = 7, .position={1.5,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});
	e.add_force_to_type(ConstantForce(3,4,5), 7);
	e.build();

	e.update_forces();
	auto const& out = e.export_particles();
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_EQ(p1.force, vec3(0,0,0));
	EXPECT_EQ(p2.force, vec3(0,0,0));
}

TEST(LinkedCellsTest, TwoParticles_IdSpecificForce) {
    Environment e;
	e.set_container(std::make_unique<LinkedCells>());

    e.add_particle(Particle{.id = 42, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 99, .type = 0, .position={0,1,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_force_to_type(NoForce(), 0);

    e.add_force_between_ids(ConstantForce(-1,2,-3), 42, 99);
    e.build();

    e.update_forces();
    auto const& out = e.export_particles();
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TEST(LinkedCellsTest, TwoParticles_InverseSquare) {
    Environment e;
	e.set_container(std::make_unique<LinkedCells>());

    e.set_extent({10,10,10});

	e.add_particle(Particle{.id = 0, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add_particle(Particle{.id = 1, .type = 1, .position={2,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});

	e.add_force_to_type(NoForce(), 0);
	e.add_force_to_type(NoForce(), 1);

    e.add_force_between_types(InverseSquare(5.0), 0, 1);
    e.build();

    e.update_forces();
    auto const& out = e.export_particles();
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


class OrbitMonitor final : public io::Monitor {
public:
	OrbitMonitor(): Monitor(1) {}
	explicit OrbitMonitor(const double v, const double r): Monitor(1), v(v), r(r) {}

	void record(size_t i, double, const std::vector<env::impl::Particle>& particles) const {
		const auto p = particles[0].mass < 1 ?  particles[0]: particles[1];

		EXPECT_NEAR(p.velocity.norm(), v, 1e-3) << "Velocity mismatch at step " << i;
		EXPECT_NEAR(p.position.norm(), r, 1e-3) << "Position mismatch at step " << i;
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

	Environment env;
	env.add_particle({0,R,0}, {v, 0, 0}, m);
	env.add_particle({0,0,0}, {0, 0, 0}, M);
	env.add_force_to_type(InverseSquare(G), 0);
	env.set_container(std::make_unique<LinkedCells>(v));
	env.set_origin({-1.5*v,-1.5*v,0});
	env.set_extent({3*v,3*v,1});

	StoermerVerlet<OrbitMonitor> integrator(env);
	integrator.add_monitor<OrbitMonitor>(OrbitMonitor(v, R));
	integrator.run(0.001, T);

	std::vector<env::impl::Particle> particles;
	for (auto & p : env.particles()) {
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


