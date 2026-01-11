#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <utils.h>

#include "april/april.hpp"

using namespace april;

// Sinks to store test data
struct SpyFieldSinks {
    int init_call_count = 0;
    int update_call_count = 0;
    int apply_call_count = 0; // Will be N_steps * N_particles
    std::vector<size_t> steps_at_update;
};

// A SpyField that uses a raw pointer to sinks owned by the test
class SpyField : public field::Field {
public:
    static constexpr env::FieldMask fields = to_field_mask(env::Field::force);
    SpyFieldSinks* sinks = nullptr;

    explicit SpyField(SpyFieldSinks* sinks_ptr) : sinks(sinks_ptr) {}
    SpyField() {} // Default ctor is required

    template<class S>
    void init(const core::SystemContext<S>&) {
        if (sinks) sinks->init_call_count++;
    }

    template<class S>
    void update(const core::SystemContext<S>& ctx) {
        if (sinks) {
            sinks->update_call_count++;
            sinks->steps_at_update.push_back(ctx.step());
        }
    }

    template<env::IsUserData U>
    void apply(const env::RestrictedParticleRef<fields, U>& /*particle*/) const {
        if (sinks) {
            sinks->apply_call_count++;
        }
    }
};


class FieldTest : public testing::Test {
protected:
    SpyFieldSinks sinks;

    // Helper to set up a minimal environment with N particles
    auto setup_environment(const int num_particles) {
        sinks = {}; // Reset sinks

        Environment env(
            forces<NoForce>,
            controllers<>,
            fields<SpyField> // Register our SpyField type
        );

        for (int i = 0; i < num_particles; ++i) {
            env.with_particle(Particle().at({static_cast<vec3::type>(i),0,0}).as_type(0).with_mass(1));
        }

        // Chain the rest of the setup
        return env.with_force(NoForce(), to_type(0))
                  .with_extent(10,10,10)
                  .with_field(SpyField(&sinks));
    }
};


TEST_F(FieldTest, InitIsCalledOnce) {
    const auto env = setup_environment(2);
    auto system = build_system(env, DirectSumAoS());
    VelocityVerlet(system).with_dt(0.01).for_steps(10).run();

    EXPECT_EQ(sinks.init_call_count, 1);
}

TEST_F(FieldTest, UpdateIsCalledEveryStep) {
    const auto env = setup_environment(2);
    auto system = build_system(env, DirectSumAoS());
    VelocityVerlet(system).with_dt(0.01).for_steps(5).run(); // Runs 0, 1, 2, 3, 4

    EXPECT_EQ(sinks.update_call_count, 5);
    std::vector<size_t> expected_steps = {0, 1, 2, 3, 4};
    EXPECT_EQ(sinks.steps_at_update, expected_steps);
}

TEST_F(FieldTest, ApplyIsCalledPerParticlePerStep) {
    const int num_particles = 3;
    const int num_steps = 5;

    const auto env = setup_environment(num_particles);
    auto system = build_system(env, DirectSumAoS());
    VelocityVerlet(system).with_dt(0.01).for_steps(num_steps).run(); // Runs 0, 1, 2, 3, 4

    EXPECT_EQ(sinks.apply_call_count, num_particles * num_steps);
    EXPECT_EQ(sinks.init_call_count, 1);
    EXPECT_EQ(sinks.update_call_count, num_steps);
}

// A second spy field to test multi-type field registration
class SpyField2 : public field::Field {
public:
    static constexpr env::FieldMask fields = to_field_mask(env::Field::force);
    SpyFieldSinks* sinks = nullptr;
    explicit SpyField2(SpyFieldSinks* sinks_ptr) : sinks(sinks_ptr) {}
    SpyField2() : sinks(nullptr) {}

    template<class S>
    void init(const core::SystemContext<S>&) {
        if (sinks) sinks->init_call_count++;
    }
    template<class S>
    void update(const core::SystemContext<S>&) {
        if (sinks) sinks->update_call_count++;
    }
    template<env::IsUserData U>
    void apply(const env::RestrictedParticleRef<fields, U>&) const {
        if (sinks) sinks->apply_call_count++;
    }
};


TEST_F(FieldTest, MultipleDifferentSpyFields) {
    SpyFieldSinks sinks1;
    SpyFieldSinks sinks2;

    const int num_particles = 2;
    const int num_steps = 5;

    Environment env(
        forces<NoForce>,
        controllers<>,
        fields<SpyField, SpyField2> // Register both types
    );

    for (int i = 0; i < num_particles; ++i) {
        env.with_particle(Particle().at({static_cast<vec3::type>(i),0,0}).as_type(0).with_mass(1));
    }

    // Use chained API
    env.with_force(NoForce(), to_type(0))
       .with_extent(10,10,10)
       .with_field(SpyField(&sinks1))
       .with_field(SpyField2(&sinks2));

    auto system = build_system(env, DirectSumAoS());
    VelocityVerlet(system).with_dt(0.01).for_steps(num_steps).run();

    // Check sinks for SpyField 1
    EXPECT_EQ(sinks1.init_call_count, 1);
    EXPECT_EQ(sinks1.update_call_count, num_steps);
    EXPECT_EQ(sinks1.apply_call_count, num_particles * num_steps);

    // Check sinks for SpyField 2
    EXPECT_EQ(sinks2.init_call_count, 1);
    EXPECT_EQ(sinks2.update_call_count, num_steps);
    EXPECT_EQ(sinks2.apply_call_count, num_particles * num_steps);
}


TEST(FieldIntegrationTest, UniformFieldModifiesForce) {
    const vec3 field_force = {5.0, 6.0, 7.0};

    // Use chained API
    Environment env(
        forces<NoForce>,
        controllers<>,
        fields<UniformField> // Register UniformField
    );
    env.with_particle(Particle().at({1,0,0}).as_type(0).with_mass(1))
       .with_particle(Particle().at({2,0,0}).as_type(0).with_mass(1))
       .with_force(NoForce(), to_type(0))
       .with_extent(10,10,10)
       .with_field(UniformField(field_force));

    auto system = build_system(env, DirectSumAoS());

    VelocityVerlet(system).with_dt(0.01).for_steps(1).run();

    auto particles = export_particles(system);
    ASSERT_EQ(particles.size(), 2);

    EXPECT_NEAR(particles[0].force.x, field_force.x, 1e-12);
    EXPECT_NEAR(particles[0].force.y, field_force.y, 1e-12);
    EXPECT_NEAR(particles[0].force.z, field_force.z, 1e-12);

    EXPECT_NEAR(particles[1].force.x, field_force.x, 1e-12);
    EXPECT_NEAR(particles[1].force.y, field_force.y, 1e-12);
    EXPECT_NEAR(particles[1].force.z, field_force.z, 1e-12);
}


TEST(FieldIntegrationTest, MultipleDifferentFieldsAreAdditive) {
    const vec3 uniform_force = {1.0, 1.0, 1.0};
    const vec3 local_force = {10.0, 0.0, 0.0};

    env::Domain local_region;
    local_region.origin = {0,0,0};
    local_region.extent = {5,5,5};

    // Use chained API
    Environment env(
        forces<NoForce>,
        controllers<>,
        fields<UniformField, LocalForceField> // Register both types
    );
    env.with_particle(Particle().at({1,1,1}).as_type(0).with_mass(1).with_id(1))
       .with_particle(Particle().at({9,9,9}).as_type(0).with_mass(1).with_id(2))
       .with_force(NoForce(), to_type(0))
       .with_extent(10,10,10)
       .with_field(UniformField(uniform_force))
       .with_field(LocalForceField(local_force, local_region, 0.0, 99.0));

    BuildInfo info;
    auto system = build_system(env, DirectSumAoS(), &info);
    VelocityVerlet(system).with_dt(0.01).for_steps(1).run(); // Runs step 0

    auto particles = export_particles(system);
    ASSERT_EQ(particles.size(), 2);

    auto p1 = (particles[0].id == info.id_map[1]) ? particles[0] : particles[1];
    auto p2 = (particles[0].id == info.id_map[1]) ? particles[1] : particles[0];

    // Particle 1 should have BOTH forces
    EXPECT_NEAR(p1.force.x, uniform_force.x + local_force.x, 1e-12);
    EXPECT_NEAR(p1.force.y, uniform_force.y + local_force.y, 1e-12);
    EXPECT_NEAR(p1.force.z, uniform_force.z + local_force.z, 1e-12);

    // Particle 2 should have ONLY the uniform force
    EXPECT_NEAR(p2.force.x, uniform_force.x, 1e-12);
    EXPECT_NEAR(p2.force.y, uniform_force.y, 1e-12);
    EXPECT_NEAR(p2.force.z, uniform_force.z, 1e-12);
}