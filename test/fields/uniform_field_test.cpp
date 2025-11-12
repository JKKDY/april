#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;

// Helper to create a dummy particle record for testing
// We don't need the full utils.h for this.
static env::internal::ParticleRecord<env::NoUserData>
make_test_particle(const vec3& force) {
	env::internal::ParticleRecord<env::NoUserData> p;
	p.force = force;
	// Other fields are not needed by UniformField,
	// but we set them for completeness.
	p.id = 0;
	p.type = 0;
	p.mass = 1.0;
	p.state = env::ParticleState::ALIVE;
	p.position = {0,0,0};
	p.velocity = {0,0,0};
	return p;
}


TEST(UniformFieldTest, ApplyIsAdditive) {
	// 1. Create the field
	const vec3 field_force = {1.0, 2.0, 3.0};
	const auto field = field::UniformField(field_force);

	// 2. Create a mock particle with an initial force
	auto p_rec = make_test_particle({10.0, 0.0, 0.0});

	// 3. Create the necessary reference for the 'apply' method
	auto p_fetcher = env::internal::ParticleRecordFetcher(p_rec);
	auto p_ref = env::RestrictedParticleRef<field::UniformField::fields, env::NoUserData>(p_fetcher);

	// 4. Call apply()
	field.apply(p_ref);

	// 5. Check that the force was added
	EXPECT_NEAR(p_rec.force.x, 10.0 + field_force.x, 1e-12);
	EXPECT_NEAR(p_rec.force.y,  0.0 + field_force.y, 1e-12);
	EXPECT_NEAR(p_rec.force.z,  0.0 + field_force.z, 1e-12);

	// 6. Call apply() again
	field.apply(p_ref);

	// 7. Check that the force was added again
	EXPECT_NEAR(p_rec.force.x, 10.0 + (2 * field_force.x), 1e-12);
	EXPECT_NEAR(p_rec.force.y,  0.0 + (2 * field_force.y), 1e-12);
	EXPECT_NEAR(p_rec.force.z,  0.0 + (2 * field_force.z), 1e-12);
}