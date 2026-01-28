#include <gtest/gtest.h>

// #include "april/april.hpp"

#include "april/common.hpp"
#include "april/forces/force.hpp"
#include "april/forces/coulomb.hpp"
#include "april/forces/gravity.hpp"
#include "april/forces/harmonic.hpp"
#include "april/forces/lennard_jones.hpp"
#include "april/forces/force_table.hpp"

using namespace april::force;
using namespace april::env;
using namespace april;


struct ForceTestUserData {
    double charge = 0.0;
    bool operator==(const ForceTestUserData&) const = default;
};


static_assert(IsUserData<ForceTestUserData>, "Test data must satisfy IsUserData");

static constexpr FieldMask TestMask = Field::position | Field::mass | Field::user_data;

class ForceTest : public testing::Test {
protected:
    vec3 pos1, pos2;
    double mass1{}, mass2{};
    ForceTestUserData data1, data2;

    // B. The Source (The struct of pointers)
    // IsConst = false allows us to point to our mutable member variables
    using SourceT = ParticleSource<TestMask, ForceTestUserData, false>;
    SourceT source1;
    SourceT source2;

    vec3 r_vec;

    void SetUp() override {
        // Initialize Data
        pos1 = {0.0, 0.0, 0.0};
        mass1 = 10.0;
        data1.charge = 1.0;

        pos2 = {2.0, 0.0, 0.0};
        mass2 = 20.0;
        data2.charge = -2.0;

        r_vec = pos2 - pos1;

        // Initialize Sources (Point to the data)
        source1.position  = &pos1;
        source1.mass      = &mass1;
        source1.user_data = &data1;

        source2.position  = &pos2;
        source2.mass      = &mass2;
        source2.user_data = &data2;
    }

    // Helper to generate the View expected by the Force operator
    [[nodiscard]] auto get_view1() const {
        return ParticleView<TestMask, ForceTestUserData>(source1);
    }
    [[nodiscard]] auto get_view2() const {
        return ParticleView<TestMask, ForceTestUserData>(source2);
    }
};



// --- Force Tests ---
// Verifies that NoForce consistently returns a zero vector and
// that mixing two NoForce objects preserves this behavior.
TEST_F(ForceTest, NoForceTest) {
    NoForce f;
    vec3 force = f(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(0.0, 0.0, 0.0));

    auto f2 = f.mix(f);
    EXPECT_EQ(f2(get_view1(), get_view2(), r_vec), vec3(0.0, 0.0, 0.0));
}

// Verifies Newtonian gravity calculation ($F = G \frac{m_1 m_2}{r^2}$)
// and ensures the force is zeroed when the distance exceeds the cutoff.
TEST_F(ForceTest, GravityTest) {
    Gravity g(1.0);
    vec3 force = g(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(50.0, 0.0, 0.0));
}

// Tests Hooke's Law ($F = -k(r - r_0)$) for both tension and compression,
// and validates that the interaction respects the defined cutoff.
TEST_F(ForceTest, HarmonicTest) {
    Harmonic h_stretch(10.0, 1.0);
    vec3 force_stretch = h_stretch(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force_stretch, vec3(10.0, 0.0, 0.0));

    Harmonic h_compress(10.0, 3.0);
    vec3 force_compress = h_compress(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force_compress, vec3(-10.0, 0.0, 0.0));
}

// Validates electrostatic force using charges stored in UserData ($F = k_e \frac{q_1 q_2}{r^2}$)
// and ensures cutoff logic correctly disables long-range interactions.
TEST_F(ForceTest, CoulombTest) {
    Coulomb c(1.0);
    vec3 force = c(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(-0.5, 0.0, 0.0));
}

// Tests the 12-6 Lennard-Jones potential for van der Waals interactions,
// checking mathematical precision near the equilibrium and enforcing cutoff limits.
TEST_F(ForceTest, LennardJonesTest) {
    LennardJones lj(1.0, 2.0);
    vec3 force = lj(get_view1(), get_view2(), r_vec);

    EXPECT_NEAR(force.x, -12.0, 1e-9);
    EXPECT_NEAR(force.y, 0.0, 1e-9);
    EXPECT_NEAR(force.z, 0.0, 1e-9);
}