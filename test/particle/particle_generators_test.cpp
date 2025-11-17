#include <gtest/gtest.h>
#include <string>
#include <any>
#include <stdexcept>

#include "april/particle/generators.h"
#include "april/particle/particle.h"
#include "april/common.h"
#include "april/particle/defs.h"

using namespace april::env;
using namespace april;


// --- ParticleCuboid Tests ---
TEST(ParticleGeneratorTest, CuboidHappyPath) {
    const std::any test_data = std::string("cuboid_test");

    const auto gen = ParticleCuboid()
        .at(1.0, 2.0, 3.0)
        .velocity(10.0, 0.0, 0.0)
        .count(2, 3, 4) // 24 particles total
        .spacing(0.5)
        .mass(1.2)
        .type(7)
        .state(ParticleState::ALIVE)
        .with_data(test_data);

    std::vector<Particle> particles = gen.to_particles();

    ASSERT_EQ(particles.size(), 24);

    // Check first particle (at origin)
    const auto& p_first = particles.front();
    EXPECT_EQ(p_first.position, vec3(1.0, 2.0, 3.0)); // origin
    EXPECT_EQ(p_first.velocity, vec3(10.0, 0.0, 0.0));
    EXPECT_EQ(p_first.mass, 1.2);
    EXPECT_EQ(p_first.type, 7);
    EXPECT_EQ(p_first.state, ParticleState::ALIVE);
    ASSERT_NO_THROW({auto _ = std::any_cast<std::string>(p_first.user_data);});
    EXPECT_EQ(std::any_cast<std::string>(p_first.user_data), "cuboid_test");

    // Check last particle (at max coordinate)
    // Count is (2,3,4) -> max index is (1,2,3)
    // Position = origin + index * spacing
    // Position = {1,2,3} + {1,2,3} * 0.5 = {1,2,3} + {0.5, 1.0, 1.5} = {1.5, 3.0, 4.5}
    const auto& p_last = particles.back();
    EXPECT_EQ(p_last.position, vec3(1.5, 3.0, 4.5));
    EXPECT_EQ(p_last.mass, 1.2);
}

TEST(ParticleGeneratorTest, CuboidOverloads) {
    // Use the vec3/uint3 overloads
    const auto gen = ParticleCuboid()
        .at(vec3(1, 2, 3))
       .velocity(vec3(4, 5, 6))
       .count(uint3{2, 2, 2})
       .spacing(1.0)
       .mass(1.0);

    const std::vector<Particle> particles = gen.to_particles();
    ASSERT_EQ(particles.size(), 8);
    EXPECT_EQ(particles.front().position, vec3(1, 2, 3));
    EXPECT_EQ(particles.front().velocity, vec3(4, 5, 6));
}

TEST(ParticleGeneratorTest, CuboidErrorZeroDistance) {
    ParticleCuboid gen;
    gen.count(1, 1, 1).spacing(0.0); // Set distance to 0


    EXPECT_THROW({const auto _ = gen.to_particles();}, std::logic_error);
}



// --- ParticleSphere Tests ---
TEST(ParticleGeneratorTest, SphereHappyPathOneParticle) {
    // Create a sphere where the spacing is larger than the radius.
    // This should result in eff_radii = spacing, and the loops
    // should only allow the single particle at (0,0,0) to pass.
    const std::any test_data = std::string("sphere_test");

    const auto gen = ParticleSphere()
        .at(10.0, 10.0, 10.0)
       .velocity(1.0, 2.0, 3.0)
       .radius(0.5) // radius < spacing
       .spacing(1.0)
       .mass(5.0)
       .type(3)
       .state(ParticleState::ALIVE)
       .with_data(test_data);
       
    std::vector<Particle> particles = gen.to_particles();

    // Only the center particle (0,0,0) should be generated
    ASSERT_EQ(particles.size(), 1);
    
    const auto& p = particles.front();
    EXPECT_EQ(p.position, vec3(10.0, 10.0, 10.0)); // center
    EXPECT_EQ(p.velocity, vec3(1.0, 2.0, 3.0));
    EXPECT_EQ(p.mass, 5.0);
    EXPECT_EQ(p.type, 3);
    EXPECT_EQ(p.state, ParticleState::ALIVE);
    ASSERT_NO_THROW({auto _ = std::any_cast<std::string>(p.user_data); });
    EXPECT_EQ(std::any_cast<std::string>(p.user_data), "sphere_test");
}

TEST(ParticleGeneratorTest, SphereOverloads) {
	// Use the vec3 overloads
	const auto gen = ParticleSphere()
        .at(vec3(1, 2, 3))
        .velocity(vec3(4, 5, 6))
        .radius_xyz(vec3(10, 10, 10)) // Use ellipsoid overload
        .spacing(100.0) // Ensure only 1 particle
        .mass(1.0);

    std::vector<Particle> particles = gen.to_particles();
    ASSERT_EQ(particles.size(), 1);
    EXPECT_EQ(particles.front().position, vec3(1, 2, 3));
    EXPECT_EQ(particles.front().velocity, vec3(4, 5, 6));
}

TEST(ParticleGeneratorTest, SphereErrorZeroDistance) {
    ParticleSphere gen;
    gen.radius(1.0).spacing(0.0); // Set distance to 0

    EXPECT_THROW({const auto _ = gen.to_particles();}, std::logic_error);
}

TEST(ParticleGeneratorTest, Sphere2DGeneration) {
    // 2D circle in the XY plane
    const auto gen = ParticleSphere()
        .at(10.0, 10.0, 10.0)
        .radius_xyz(2.0, 2.0, 0.0) // radii.z = 0
        .spacing(1.0)
        .mass(1.0);

    std::vector<Particle> particles = gen.to_particles();


    // The loop logic will check: x*x + y*y < (eff_radii.x * eff_radii.x)
    // eff_radii.x = 2.0, spacing = 1.0. Loop for x/y is from -2 to 1.
    // We check x*x + y*y < 4
    //
    // y=-2: x=-2(F), x=-1(F), x=0(F), x=1(F)
    // y=-1: x=-2(F), x=-1(T), x=0(T), x=1(T) -> 3 particles
    // y= 0: x=-2(F), x=-1(T), x=0(T), x=1(T) -> 3 particles
    // y= 1: x=-2(F), x=-1(T), x=0(T), x=1(T) -> 3 particles
    // Total = 9 particles

    ASSERT_EQ(particles.size(), 9);

    // Also check that all particles are in the Z=0 plane
    const double expected_z = gen.center.z;
    for (const auto& p : particles) {
        EXPECT_EQ(p.position.z, expected_z);
    }
}



// --- Common Feature Tests ---
TEST(ParticleGeneratorTest, ThermalVelocity) {
    // thermal function: v_thermal = position
    auto thermal_fn = [](const vec3& pos) { return pos; };

    auto gen = ParticleCuboid()
        .at(5, 5, 5)
        .velocity(1, 1, 1) // mean velocity
        .count(1, 1, 1)
        .spacing(1.0)
        .mass(1.0)
        .thermal(thermal_fn);

    std::vector<Particle> particles = gen.to_particles();

    ASSERT_EQ(particles.size(), 1);
    
    // check default thermal (should be zero)
    ParticleCuboid gen_default;
    gen_default.count(1,1,1).velocity(1,1,1).spacing(1.0);
    auto p_default = gen_default.to_particles().front();
    EXPECT_EQ(p_default.velocity, vec3(1,1,1)); // velocity = mean + 0

    // check custom thermal
    // final velocity = mean_velocity + thermal_fn(position)
    // final velocity = {1,1,1} + {5,5,5} = {6,6,6}
    const auto& p_thermal = particles.front();
    EXPECT_EQ(p_thermal.position, vec3(5, 5, 5));
    EXPECT_EQ(p_thermal.velocity, vec3(6, 6, 6));
}