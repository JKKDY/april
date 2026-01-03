#include <utils.h>
#include <gtest/gtest.h>
#include "april/april.hpp"

using namespace april;

TEST(LocalFieldTest, SpatialCheck) {
    // Test that the field applies force ONLY
    // to particles inside its defined region.
    
    const vec3 field_force = {10.0, 0.0, 0.0};
    
    // Define a local region from {5,5,5} to {10,10,10}
    env::Domain local_region;
    local_region.origin = {5,5,5};
    local_region.extent = {5,5,5};

    Environment env(
        forces<NoForce>,
        controllers<>,
        fields<LocalForceField>
    );
    
    // p1 is inside the region
    env.with_particle(Particle().at({7, 7, 7}).as_type(0).with_mass(1).with_id(1));
    // p2 is outside the region
    env.with_particle(Particle().at({1, 1, 1}).as_type(0).with_mass(1).with_id(2));
    
    env.with_force(NoForce(), to_type(0))
       .with_extent(20, 20, 20)
       // Add the field, active from t=0.0 to t=10.0
       .with_field(LocalForceField(field_force, local_region, 0.0, 10.0));

    BuildInfo info;
    auto system = build_system(env, DirectSumAoS(), &info);
    system.apply_force_fields();
    // Run at step 0 (t=0.0), which is inside the active time
    // StoermerVerlet(system).with_dt(0.01).for_steps(1).run();

    auto particles = export_particles(system);
    ASSERT_EQ(particles.size(), 2);

    auto p1 = (particles[0].id == info.id_map[1]) ? particles[0] : particles[1];
    auto p2 = (particles[0].id == info.id_map[1]) ? particles[1] : particles[0];

    // p1 (inside) should have the force
    EXPECT_NEAR(p1.force.x, field_force.x, 1e-12);
    
    // p2 (outside) should have no force
    EXPECT_NEAR(p2.force.x, 0.0, 1e-12);
}


TEST(LocalFieldTest, TimeCheck) {
    // Test that the field turns on and off
    // at the correct simulation times.
    
    const vec3 field_force = {10.0, 0.0, 0.0};
    
    // Region covers the whole domain
    env::Domain local_region;
    local_region.origin = {0,0,0};
    local_region.extent = {20,20,20};

    Environment env(
        forces<NoForce>,
        controllers<>,
        fields<LocalForceField>
    );
    
    // Particle is always inside the region
    env.with_particle(Particle().at({5, 5, 5}).as_type(0).with_mass(1).with_id(1));
    
    env.with_force(NoForce(), to_type(0))
       .with_extent(20, 20, 20)
       // Field is active ONLY between t=0.025 and t=0.045
       .with_field(LocalForceField(field_force, local_region, 0.025, 0.045));

    auto system = build_system(env, DirectSumAoS());
    auto integrator = VelocityVerlet(system).with_dt(0.01);

    // Phase 1: before (t=0.0, t=0.01)
    integrator.for_steps(2).run(); // Runs step 0 (t=0.0), 1 (t=0.01)
    auto p1 = export_particles(system)[0];
    EXPECT_NEAR(p1.force.x, 0.0, 1e-12); // Should not be active yet

    // Phase 2: during (t=0.02, t=0.03)
    integrator.for_steps(2).run(); // Runs step 2 (t=0.02), 3 (t=0.03)
    auto p2 = export_particles(system)[0];
    // Step 2 (t=0.02): update() -> active=false.
    // Step 3 (t=0.03): update() -> active=true.
    EXPECT_NEAR(p2.force.x, field_force.x, 1e-12); // Should be active

    // Phase 3: after (t=0.04, t=0.05)
    integrator.for_steps(2).run(); // Runs step 4 (t=0.04), 5 (t=0.05)
    auto p3 = export_particles(system)[0];
    // Step 4 (t=0.04): update() -> active=true.
    // Step 5 (t=0.05): update() -> active=false.
    EXPECT_NEAR(p3.force.x, 0.0, 1e-12); // Should be inactive again
}