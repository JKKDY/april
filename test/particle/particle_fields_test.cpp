#include <gtest/gtest.h>
#include <string>
#include <type_traits> // std::is_same_v

#include "april/base/types.hpp"
#include "april/particle/scalar_access.hpp"
#include "april/particle/particle.hpp"

using namespace april::particle::internal;
using namespace april;

struct MyTestUserData {
    int id = 0;
    double value = 0.0;

    // Add an equality operator for GTest assertions
    bool operator==(const MyTestUserData&) const = default;
};

using TestUserDataT = MyTestUserData;

// make sure test type is valid
static_assert(core::IsParticleAttributes<TestUserDataT>, "MyTestUserData does not satisfy IsUserData");


// create a particle record with data
class ParticleViewsTest : public testing::Test {
protected:
    // The backing storage
    particle::ParticleRecord<TestUserDataT> particle_data;

    void SetUp() override {
        particle_data.id = 123;
        particle_data.type = 4;
        particle_data.position = {1.0, 2.0, 3.0};
        particle_data.velocity = {4.0, 5.0, 6.0};
        particle_data.force = {7.0, 8.0, 9.0};
        particle_data.old_position = {10.0, 11.0, 12.0};
        particle_data.mass = 1.1;
        particle_data.state = ParticleState::ALIVE;
        particle_data.attributes = MyTestUserData{10, 20.5};
    }

    // Helper: Create a Mutable Source pointing to the Record's fields
    auto get_source() {
        particle::internal::ParticleSource<ParticleField::all, TestUserDataT, false> src;
        src.position     = &particle_data.position;
        src.velocity     = &particle_data.velocity;
        src.force        = &particle_data.force;
        src.old_position = &particle_data.old_position;
        src.mass         = &particle_data.mass;
        src.state        = &particle_data.state;
        src.type         = &particle_data.type;
        src.id           = &particle_data.id;
        src.attributes    = &particle_data.attributes;
        return src;
    }

    // Helper: Create a Const Source pointing to the Record's fields
    auto get_const_source() {
        particle::internal::ParticleSource<ParticleField::all, TestUserDataT, true> src;
        src.position     = &particle_data.position;
        src.velocity     = &particle_data.velocity;
        src.force        = &particle_data.force;
        src.old_position = &particle_data.old_position;
        src.mass         = &particle_data.mass;
        src.state        = &particle_data.state;
        src.type         = &particle_data.type;
        src.id           = &particle_data.id;
        src.attributes    = &particle_data.attributes;
        return src;
    }
};


TEST(ParticleViewsHelpersTest, BitmaskOperators) {
    constexpr auto mask1 = ParticleField::position | ParticleField::velocity;
    EXPECT_EQ(static_cast<size_t>(mask1), (1u << 0) | (1u << 1));

    constexpr ParticleField mask2 = mask1 | ParticleField::force;
    EXPECT_EQ(static_cast<size_t>(mask2), (1u << 0) | (1u << 1) | (1u << 2));

    constexpr ParticleField mask3 = ParticleField::id | mask2;
    EXPECT_EQ(static_cast<size_t>(mask3), (1u << 0) | (1u << 1) | (1u << 2) | (1u << 7));

    // Test has_field_v
    EXPECT_TRUE(static_cast<bool>(mask3 & ParticleField::position));
    EXPECT_TRUE(static_cast<bool>(mask3 & ParticleField::id));
    EXPECT_FALSE(static_cast<bool>(mask3 & ParticleField::mass));
    EXPECT_TRUE(static_cast<bool>(ParticleField::all & ParticleField::attributes));
    EXPECT_FALSE(static_cast<bool>(ParticleField::none & ParticleField::position));
}


// --- Test ParticleRef ---
TEST_F(ParticleViewsTest, ParticleRefAllFieldsRead) {
    auto src = get_source();
    const ScalarParticleRef<ParticleField::all, TestUserDataT> ref(src);

    EXPECT_EQ(ref.position, particle_data.position);
    EXPECT_EQ(ref.velocity, particle_data.velocity);
    EXPECT_EQ(ref.force, particle_data.force);
    EXPECT_EQ(ref.old_position, particle_data.old_position);
    EXPECT_EQ(ref.mass, particle_data.mass);
    EXPECT_EQ(ref.state, particle_data.state);
    EXPECT_EQ(ref.type, particle_data.type); // Read-only copy
    EXPECT_EQ(ref.id, particle_data.id);     // Read-only copy
    EXPECT_EQ(ref.attributes, particle_data.attributes);
}

TEST_F(ParticleViewsTest, ParticleRefAllFieldsWrite) {
    const auto src = get_source();
    ScalarParticleRef<ParticleField::all, TestUserDataT> ref(src);

    constexpr MyTestUserData updated_data{99, -1.0};

    // modify data
    ref.position = {101.0, 102.0, 103.0};
    ref.mass = 2.2;
    ref.attributes = updated_data;

    // check for modification
    EXPECT_EQ(particle_data.position, vec3(101.0, 102.0, 103.0));
    EXPECT_DOUBLE_EQ(particle_data.mass, 2.2);
    EXPECT_EQ(particle_data.attributes, updated_data);
}

TEST_F(ParticleViewsTest, ParticleRefPartialMask) {
    constexpr auto mask = ParticleField::position | ParticleField::mass | ParticleField::attributes;

    auto src = get_source(); // Source has ALL fields populated
    ScalarParticleRef<mask, TestUserDataT> ref(src); // Ref only maps subset

    // check present fields are correct
    EXPECT_EQ(ref.position, particle_data.position);
    EXPECT_TRUE((std::is_same_v<decltype(ref.position), math::Vec3Proxy<vec3::type>>));
    EXPECT_EQ(ref.mass, particle_data.mass);
    EXPECT_TRUE((std::is_same_v<decltype(ref.mass), double&>));
    EXPECT_EQ(ref.attributes, particle_data.attributes);
    EXPECT_TRUE((std::is_same_v<decltype(ref.attributes), TestUserDataT&>));

    // check that absent fields are absent
    EXPECT_TRUE(!(std::is_same_v<decltype(ref.velocity), math::Vec3Proxy<vec3::type>>));
    EXPECT_TRUE(!(std::is_same_v<decltype(ref.force), math::Vec3Proxy<vec3::type>>));
}


// --- Test ParticleView ---
TEST_F(ParticleViewsTest, ParticleViewIsConst) {
    auto src = get_const_source(); // Source is const
    ScalarParticleView<ParticleField::all, TestUserDataT> view(src);

    // check values
    EXPECT_EQ(view.position, particle_data.position);
    EXPECT_EQ(view.mass, particle_data.mass);
    EXPECT_EQ(view.attributes, particle_data.attributes);

    // check types are const (or copies)
    EXPECT_TRUE((std::is_same_v<decltype(view.position), const math::Vec3Proxy<const vec3::type>>));
    EXPECT_TRUE((std::is_same_v<decltype(view.mass), const double&>)); // copy
    EXPECT_TRUE((std::is_same_v<decltype(view.attributes), const TestUserDataT&>));
}


// --- Test RestrictedParticleRef ---
TEST_F(ParticleViewsTest, RestrictedParticleRefAccess) {
    constexpr auto mask = ParticleField::position | ParticleField::force | ParticleField::id | ParticleField::attributes;
    auto src = get_source(); // Mutable source

    ScalarRestrictedParticleRef<mask, TestUserDataT> restricted_ref(src);

    // check that 'force' is mutable
    EXPECT_TRUE((std::is_same_v<decltype(restricted_ref.force), math::Vec3Proxy<vec3::type>>));

    // check that other fields are const or copies
    EXPECT_TRUE((std::is_same_v<decltype(restricted_ref.position), const math::Vec3Proxy<const vec3::type>>));
    EXPECT_TRUE((std::is_same_v<decltype(restricted_ref.id), const ParticleID>)); // copy
    EXPECT_TRUE((std::is_same_v<decltype(restricted_ref.attributes), const TestUserDataT&>));

    // check that absent fields are monostate
    EXPECT_TRUE((!std::is_same_v<decltype(restricted_ref.velocity), math::Vec3Proxy<vec3::type>>));
    EXPECT_TRUE((!std::is_same_v<decltype(restricted_ref.mass), math::Vec3Proxy<vec3::type>>));

    // Test write access
    restricted_ref.force = {999.0, 999.0, 999.0};
    EXPECT_EQ(particle_data.force, vec3(999.0, 999.0, 999.0));
}





