#include <gtest/gtest.h>
#include <gmock/gmock.h>

using testing::AnyOf;
using testing::Eq;

#include "april/april.h"
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


TEST(DirectSumTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
    e.add(Particle{.id = 0, .type = 0, .position={1,2,3},.velocity={0,0,0}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));

	auto sys = core::build_system(e, DirectSum());
    sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(DirectSumTest, TwoParticles_ConstantTypeForce) {
    Environment e (forces<ConstantForce>);
	e.add(Particle{.id = 0, .type = 7, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 1, .type = 7, .position={1,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = core::build_system(e, DirectSum());
    sys.update_forces();
    auto const& out = sys.export_particles();

    ASSERT_EQ(out.size(), 2u);
    // both should see the same force vector
    EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5)))
	);

}

TEST(DirectSumTest, TwoParticles_IdSpecificForce) {
    Environment e (forces<ConstantForce, NoForce>);
    e.add(Particle{.id = 42, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 99, .type = 0, .position={0,1,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));

	auto sys = core::build_system(e, DirectSum());
	sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 2u);

	EXPECT_EQ(out[0].force, -out[1].force);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TEST(DirectSumTest, TwoParticles_InverseSquare) {
	Environment e (forces<InverseSquare, NoForce>);

	e.set_extent({10,10,10});
	e.add(Particle{.id = 0, .type = 0, .position={0,0,0},.velocity={}, .mass=1, .state=ParticleState::ALIVE});
    e.add(Particle{.id = 1, .type = 1, .position={2,0,0},.velocity={}, .mass=2, .state=ParticleState::ALIVE});

	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));
	e.add_force(InverseSquare(5.0), between_types(0, 1));

	auto sys = core::build_system(e, DirectSum());
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



