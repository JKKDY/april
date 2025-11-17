#include <gtest/gtest.h>
#include <string>
#include <any>

#include "april/particle/particle.h"
#include "april/common.h"  // vec3
#include "april/particle/defs.h"
#include "april/particle/generators.h"

using namespace april::env;
using namespace april;

TEST(ParticleTest, FluentSettersAndChaining) {
    // Define test data
    constexpr ParticleID test_id = 123;
    constexpr ParticleType test_type = 4;
    const vec3 test_pos = {1.0, 2.0, 3.0};
    const vec3 test_vel = {4.0, 5.0, 6.0};
    constexpr double test_mass = 7.0;
    constexpr auto test_state = ParticleState::ALIVE;
    const vec3 test_old_pos = {8.0, 9.0, 10.0};
    const vec3 test_old_force = {11.0, 12.0, 13.0};
    const vec3 test_force = {14.0, 15.0, 16.0};
    const std::any test_data = std::string("hello");

    // Chain all the setters together
    auto p = Particle().with_id(test_id)
     .as_type(test_type)
     .at(test_pos)
     .with_velocity(test_vel)
     .with_mass(test_mass)
     .with_state(test_state)
     .with_old_position(test_old_pos)
     .with_old_force(test_old_force)
     .with_force(test_force)
     .with_data(test_data);

	// check correctness
    ASSERT_TRUE(p.id.has_value());
    EXPECT_EQ(p.id.value(), test_id);

    EXPECT_EQ(p.type, test_type);

    EXPECT_EQ(p.position, test_pos);
    EXPECT_EQ(p.velocity, test_vel);

    EXPECT_DOUBLE_EQ(p.mass, test_mass);
    EXPECT_EQ(p.state, test_state);

    ASSERT_TRUE(p.old_position.has_value());
    EXPECT_EQ(p.old_position.value(), test_old_pos);

    ASSERT_TRUE(p.old_force.has_value());
    EXPECT_EQ(p.old_force.value(), test_old_force);

    ASSERT_TRUE(p.force.has_value());
    EXPECT_EQ(p.force.value(), test_force);

    ASSERT_NO_THROW({auto _ = std::any_cast<std::string>(p.user_data);});
    EXPECT_EQ(std::any_cast<std::string>(p.user_data), "hello");
}


TEST(ParticleTest, SetterOverloads) {
    // 1. ARRANGE
    Particle p;
    const vec3 expected_pos = {1.5, 2.5, 3.5};
    const vec3 expected_vel = {4.5, 5.5, 6.5};

    p = p.at(1.5, 2.5, 3.5);
    p = p.with_velocity(4.5, 5.5, 6.5);

    // 3. ASSERT
    EXPECT_EQ(p.position, expected_pos);
    EXPECT_EQ(p.velocity, expected_vel);
}


