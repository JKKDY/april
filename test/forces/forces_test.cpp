#include <gtest/gtest.h>

#include "april/forces/noforce.h"
#include "april/forces/gravity.h"
#include "april/forces/harmonic.h"
#include "april/forces/coulomb.h"
#include "april/forces/lennardjones.h"

#include "april/common.h"
#include "april/particle/fields.h"
#include "april/particle/defs.h"

using namespace april::force;
using namespace april::env;
using namespace april;


struct ForceTestUserData {
    double charge = 0.0;
    bool operator==(const ForceTestUserData&) const = default;
};


static_assert(IsUserData<ForceTestUserData>, "Test data must satisfy IsUserData");

struct MockConstFetcher {
    using user_data_t = ForceTestUserData;
    
    vec3 _pos;
    vec3 _vel;
    vec3 _force;
    vec3 _old_pos;
    vec3 _old_force;
    double _mass = 0.0;
    ParticleState _state {};
    ParticleType _type = 0;
    ParticleID _id = 0;
    user_data_t _data {};

    [[nodiscard]] const vec3& position() const { return _pos; }
    [[nodiscard]] const vec3& velocity() const { return _vel; }
    [[nodiscard]] const vec3& force() const { return _force; }
    [[nodiscard]] const vec3& old_position() const { return _old_pos; }
    [[nodiscard]] const vec3& old_force() const { return _old_force; }
    [[nodiscard]] double mass() const { return _mass; }
    [[nodiscard]] ParticleState state() const { return _state; }
    [[nodiscard]] ParticleType type() const { return _type; }
    [[nodiscard]] ParticleID id() const { return _id; }
    [[nodiscard]] const user_data_t& user_data() const { return _data; }
};
static_assert(IsConstFetcher<MockConstFetcher>, "Mock fetcher does not satisfy IsConstFetcher");

class ForceTest : public testing::Test {
protected:
    MockConstFetcher p1;
    MockConstFetcher p2;
    vec3 r_vec; 

    void SetUp() override {
        p1._pos = {0.0, 0.0, 0.0};
        p1._mass = 10.0;
        p1._data.charge = 1.0;

        p2._pos = {2.0, 0.0, 0.0};
        p2._mass = 20.0;
        p2._data.charge = -2.0;
        
        r_vec = p2._pos - p1._pos;
    }
};

// --- Force Tests ---

TEST_F(ForceTest, NoForceTest) {
    NoForce f;
    vec3 force = f(p1, p2, r_vec);
    EXPECT_EQ(force, vec3(0.0, 0.0, 0.0));
    
    auto f2 = f.mix(f);
    EXPECT_EQ(f2(p1, p2, r_vec), vec3(0.0, 0.0, 0.0));
}

TEST_F(ForceTest, GravityTest) {
    Gravity g(1.0);
    vec3 force = g(p1, p2, r_vec);
    EXPECT_EQ(force, vec3(50.0, 0.0, 0.0));

    Gravity g_cut(1.0, 1.0);
    vec3 force_cut = g_cut(p1, p2, r_vec);
    EXPECT_EQ(force_cut, vec3(0.0, 0.0, 0.0));
}

TEST_F(ForceTest, HarmonicTest) {
    Harmonic h_stretch(10.0, 1.0);
    vec3 force_stretch = h_stretch(p1, p2, r_vec);
    EXPECT_EQ(force_stretch, vec3(-10.0, 0.0, 0.0));

    Harmonic h_compress(10.0, 3.0);
    vec3 force_compress = h_compress(p1, p2, r_vec);
    EXPECT_EQ(force_compress, vec3(10.0, 0.0, 0.0));

    Harmonic h_cut(10.0, 1.0, 1.0);
    vec3 force_cut = h_cut(p1, p2, r_vec);
    EXPECT_EQ(force_cut, vec3(0.0, 0.0, 0.0));
}

TEST_F(ForceTest, CoulombTest) {
    Coulomb c(1.0);
    vec3 force = c(p1, p2, r_vec);
    EXPECT_EQ(force, vec3(-0.5, 0.0, 0.0));
    
    Coulomb c_cut(1.0, 1.0);
    vec3 force_cut = c_cut(p1, p2, r_vec);
    EXPECT_EQ(force_cut, vec3(0.0, 0.0, 0.0));
}

TEST_F(ForceTest, LennardJonesTest) {
    LennardJones lj(1.0, 2.0);
    vec3 force = lj(p1, p2, r_vec);
    EXPECT_NEAR(force.x, -12.0, 1e-9);
    EXPECT_NEAR(force.y, 0.0, 1e-9);
    EXPECT_NEAR(force.z, 0.0, 1e-9);

    LennardJones lj_cut(1.0, 2.0, 1.0);
    vec3 force_cut = lj_cut(p1, p2, r_vec);
    EXPECT_EQ(force_cut, vec3(0.0, 0.0, 0.0));
}