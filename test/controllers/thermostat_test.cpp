#include <utils.h>
#include <gtest/gtest.h>
#include "april/april.hpp"
#include <vector>

using namespace april;


// Helper functions to analyze the results of a simulation
using ParticleRec = env::internal::ParticleRecord<env::NoUserData>;

uint8_t get_dimensions(const env::Box& box) {
     return 3 -
        (box.extent.x == 0) -
        (box.extent.y == 0) -
        (box.extent.z == 0);
}

vec3 get_system_avg_v(const std::vector<ParticleRec>& particles) {
    if (particles.empty()) return {0,0,0};
    vec3 sum{0, 0, 0};
    for (const auto& p : particles) {
        sum += p.velocity;
    }
    return sum / static_cast<double>(particles.size());
}

double get_system_temp(const std::vector<ParticleRec>& particles, const vec3& avg_v, const env::Box& box) {
    double kinetic = 0.0;
    for (const auto& p : particles) {
        const vec3 dv = p.velocity - avg_v;
        kinetic += p.mass * dv.norm_squared();
    }

    const size_t dof = get_dimensions(box) * particles.size();
    if (dof == 0) return 0.0;
    return kinetic / static_cast<double>(dof);
}


// --- Temperature Calculation Tests ---

// here we check that we can use our above helpers to validate our system
TEST(ThermostatCalculationTest, InitialTemperatureTest1) {
    // 4 particles, avg_v = 0
    // K_thermal = 0.5*1*1^2 * 4 = 2.0
    // Total K = 2.0
    // N=4, D=2
    // dof = N*D = 8
    // T = K_total / dof * 2 = 4.0 / 8 = 0.5
    // Note: The original test expected 0.5. This implies a
    // different definition (e.g., T = K_per_particle / (D/2)).
    // We will test against our formula, which gives 0.25.
    Environment env(forces<NoForce>);
    env.with_particle(Particle().at({30, 10, 0}).with_velocity({-1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({70, 10, 0}).with_velocity({ 1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({30, 90, 0}).with_velocity({-1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({70, 90, 0}).with_velocity({ 1, 0, 0}).with_mass(1).as_type(0))
       .with_force(NoForce(), to_type(0))
       .with_extent(100, 100, 0); // 2D system

    auto system = build_system(env, DirectSumAoS());
    const auto particles = export_particles(system);
    const auto avg_v = get_system_avg_v(particles);
    const auto temp = get_system_temp(particles, avg_v, system.box());

    EXPECT_NEAR(avg_v.x, 0.0, 1e-12);
    EXPECT_NEAR(avg_v.y, 0.0, 1e-12);
    EXPECT_NEAR(avg_v.z, 0.0, 1e-12);
    EXPECT_NEAR(temp, 0.5, 1e-12); // K=2, dof=8. T = 2/8 = 0.25
}

TEST(ThermostatCalculationTest, InitialTemperatureTest2) {
    // 4 particles, avg_v = {1,0,0}
    // Thermal velocity (v - avg_v) is 0 for all.
    // Total kinetic energy (thermal) = 0
    // T = 0 / 8 = 0.0
    Environment env(forces<NoForce>);
    env.with_particle(Particle().at({30, 10, 0}).with_velocity({1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({70, 10, 0}).with_velocity({1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({30, 90, 0}).with_velocity({1, 0, 0}).with_mass(1).as_type(0))
       .with_particle(Particle().at({70, 90, 0}).with_velocity({1, 0, 0}).with_mass(1).as_type(0))
       .with_force(NoForce(), to_type(0))
       .with_extent(100, 100, 0); // 2D system

    auto system = build_system(env, DirectSumAoS());
    const auto particles = export_particles(system);
    const auto avg_v = get_system_avg_v(particles);
    const auto temp = get_system_temp(particles, avg_v, system.box());

    EXPECT_NEAR(avg_v.x, 1.0, 1e-12);
    EXPECT_NEAR(avg_v.y, 0.0, 1e-12);
    EXPECT_NEAR(avg_v.z, 0.0, 1e-12);
    EXPECT_NEAR(temp, 0.0, 1e-12);
}


// --- Thermostat Behavior Tests (Integration) ---
TEST(ThermostatBehaviorTest, SetInitialTemperature) {
    for (int t = 0; t < 100; t+= 10) {
        const auto cuboid = ParticleCuboid()
                            .at(0,0,0)
                            .spacing(1)
                            .count(10,10,10)
                            .type(0)
                            .mass(1);

        const auto env = Environment (forces<NoForce>, controllers<VelocityScalingThermostat>)
            .with_particles(cuboid)
            .with_force(NoForce(), to_type(0))
            .with_extent(100, 100, 100)
            .with_controller(VelocityScalingThermostat(t, 0, 30, Trigger::always()));

        auto system = build_system(env, DirectSumAoS());

        const auto particles = export_particles(system);
        const auto avg_v = get_system_avg_v(particles);
        const auto temp = get_system_temp(particles, avg_v, system.box());

        EXPECT_NEAR(avg_v.x, 0, 1);
        EXPECT_NEAR(avg_v.y, 0, 1);
        EXPECT_NEAR(avg_v.z, 0, 1);
        EXPECT_NEAR(temp, t, t * 0.05); // 5% deviation allowed
    }
}

TEST(ThermostatBehaviorTest, HoldingTemperature) {
    // Tests that velocities are constant if T_target == T_current

    const auto cuboid = ParticleCuboid()
                            .at(0,0,0)
                            .spacing(1)
                            .count(10,10,1)
                            .type(0)
                            .mass(1);

    const auto env = Environment (forces<NoForce>, boundaries<Reflective>, controllers<VelocityScalingThermostat>)
        .with_particles(cuboid)
        .with_boundaries(Reflective(), all_faces)
        .with_force(NoForce(), to_type(0))
        .with_extent(100, 100, 100)
        .with_controller(VelocityScalingThermostat(
             20, 20, 0.5, Trigger::every(10)
        ));

    auto system = build_system(env, DirectSumAoS());

    // Run for a few steps
    VelocityVerlet(system).with_dt(0.001).for_steps(100).run();

    const auto particles = export_particles(system);
    const auto avg_v = get_system_avg_v(particles);
    const auto temp = get_system_temp(particles, avg_v, system.box());

    EXPECT_NEAR(temp, 20.0, 1.0); // 5% tolerance
}


TEST(ThermostatBehaviorTest, CoolingSystem) {
    const auto cuboid = ParticleCuboid()
                          .at(0,0,0)
                          .spacing(1)
                          .count(10,10,1)
                          .type(0)
                          .mass(1);

    const auto env = Environment (forces<NoForce>, boundaries<Reflective>, controllers<VelocityScalingThermostat>)
        .with_particles(cuboid)
        .with_boundaries(Reflective(), all_faces)
        .with_force(NoForce(), to_type(0))
        .with_extent(100, 100, 100)
        .with_controller(VelocityScalingThermostat(
             20, 5, 10, Trigger::every(10)
         ));

    auto system = build_system(env, DirectSumAoS());

    // Run for a few steps
    VelocityVerlet(system).with_dt(0.001).for_steps(100).run();

    const auto particles = export_particles(system);
    const auto avg_v = get_system_avg_v(particles);
    const auto temp = get_system_temp(particles, avg_v, system.box());

    EXPECT_NEAR(temp, 5, 0.25); // 5% tolerance
}

TEST(ThermostatBehaviorTest, HeatingSystem) {
    const auto cuboid = ParticleCuboid()
                          .at(0,0,0)
                          .spacing(1)
                          .count(10,10,1)
                          .type(0)
                          .mass(1);

    const auto env = Environment (forces<NoForce>, boundaries<Reflective>, controllers<VelocityScalingThermostat>)
        .with_particles(cuboid)
        .with_boundaries(Reflective(), all_faces)
        .with_force(NoForce(), to_type(0))
        .with_extent(100, 100, 100)
        .with_controller(VelocityScalingThermostat(
            20, 80, 10, Trigger::every(10)
        ));

    auto system = build_system(env, DirectSumAoS());

    // Run for a few steps
    VelocityVerlet(system).with_dt(0.001).for_steps(100).run();

    const auto particles = export_particles(system);
    const auto avg_v = get_system_avg_v(particles);
    const auto temp = get_system_temp(particles, avg_v, system.box());

    EXPECT_NEAR(temp, 80, 80 * 0.05); // 5% tolerance
}



TEST(ThermostatBehaviorTest, Apply_HeatsThenCoolsWithTriggers) {
    constexpr double T_heat = 40.0;
    constexpr double T_cool = 5;

    const auto cuboid = ParticleCuboid()
                           .at(0,0,0)
                           .spacing(1)
                           .count(10,10,1)
                           .type(0)
                           .mass(1);

    auto env = Environment (forces<NoForce>, boundaries<Reflective>, controllers<VelocityScalingThermostat>)
        .with_particles(cuboid)
        .with_boundaries(Reflective(), all_faces)
        .with_force(NoForce(), to_type(0))
        .with_extent(100, 100, 100);

    // Controller 1: Heats to 40.0 between steps 0 and 19
    env.add_controller(VelocityScalingThermostat(
        controller::temperature_not_set, T_heat, 5.0, Trigger::between(0, 20)
    ));
    // Controller 2: Cools to 2.5 after step 20
    env.add_controller(VelocityScalingThermostat(
        controller::temperature_not_set, T_cool, 5.0, Trigger::after(20)
    ));

    auto system = build_system(env, DirectSumAoS());

    // Run heating phase
    VelocityVerlet(system).with_dt(0.01).for_steps(20).run();
    auto p_heated = export_particles(system);
    auto v_heated = get_system_avg_v(p_heated);
    auto T_heated = get_system_temp(p_heated, v_heated, system.box());
    EXPECT_NEAR(T_heated, T_heat, 0.1);

    // Run cooling phase
    VelocityVerlet(system).with_dt(0.01).for_steps(20).run();
    auto p_cooled = export_particles(system);
    auto v_cooled = get_system_avg_v(p_cooled);
    auto T_cooled = get_system_temp(p_cooled, v_cooled, system.box());
    EXPECT_NEAR(T_cooled, T_cool, 0.1);

}