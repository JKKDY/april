#include <gtest/gtest.h>

#include "utils.h"
#include "april/april.hpp"

using namespace april;

// Helper to create a dummy particle record for testing
// We don't need the full utils.h for this.
static particle::ParticleRecord<NoParticleAttributes>
make_test_particle(const vec3& force) {
	particle::ParticleRecord<NoParticleAttributes> p;
	p.force = force;
	// Other fields are not needed by UniformField,
	// but we set them for completeness.
	p.id = 0;
	p.type = 0;
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	p.position = {0,0,0};
	p.velocity = {0,0,0};
	return p;
}


TEST(UniformFieldTest, ApplyIsAdditive) {
	const vec3 field_force = {1.0, 2.0, 3.0};
	const auto field = UniformField(field_force);

	auto p_rec = make_test_particle({10.0, 0.0, 0.0});

	constexpr ParticleField Mask = UniformField::fields;

	auto get_field = [&]<ParticleField F>() {
		if constexpr (F == ParticleField::force)
			return &p_rec.force;
	};

	auto src = particle::internal::make_particle_source<Mask, Mask>(get_field);

	particle::internal::ScalarParticleRef<
		Mask,
		Mask,
		NoParticleAttributes
	> p_ref(src);

	field.apply(p_ref);

	EXPECT_NEAR(p_rec.force.x, 10.0 + field_force.x, 1e-12);
	EXPECT_NEAR(p_rec.force.y,  0.0 + field_force.y, 1e-12);
	EXPECT_NEAR(p_rec.force.z,  0.0 + field_force.z, 1e-12);

	field.apply(p_ref);

	EXPECT_NEAR(p_rec.force.x, 10.0 + (2 * field_force.x), 1e-12);
	EXPECT_NEAR(p_rec.force.y,  0.0 + (2 * field_force.y), 1e-12);
	EXPECT_NEAR(p_rec.force.z,  0.0 + (2 * field_force.z), 1e-12);
}











