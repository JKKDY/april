#include <gtest/gtest.h>
#include <april/env/environment.h>
#include "april/common.h"
#include "april/containers/direct_sum.h"
#include "april/env/particle.h"
#include <gmock/gmock.h>

#include "april/core/system.h"
using testing::AnyOf;
using testing::Eq;

using namespace april;
using namespace april::env;
using namespace april::cont;

struct ConstantForce final : Force {
	vec3 v;
	ConstantForce(const double x, const double y, const double z, const double cutoff = -1) : v{x,y,z} {
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


TEST(DirectSumTest, SingleParticle_NoForce) {
    Environment e;
    e.add(Particle{.id = 0, .type = 0, .position={1,2,3},.velocity={0,0,0}, .mass=1.0, .state=ParticleState::ALIVE});
	e.add_force(NoForce(), to_type(0));

	auto sys = core::build_system(e, DirectSum());
    sys.update_forces();

    auto const& out = sys.export_particles();
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TEST(DirectSumTest, TwoParticles_ConstantTypeForce) {
    Environment e;
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
    Environment e;
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
    Environment e;

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



