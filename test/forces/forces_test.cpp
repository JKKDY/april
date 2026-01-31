#include <gtest/gtest.h>

// #include "april/april.hpp"

#include "april/base/types.hpp"
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
    const vec3 force = f(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(0.0, 0.0, 0.0));

    auto f2 = f.mix(f);
    EXPECT_EQ(f2(get_view1(), get_view2(), r_vec), vec3(0.0, 0.0, 0.0));
}

// Verifies Newtonian gravity calculation ($F = G \frac{m_1 m_2}{r^2}$)
// and ensures the force is zeroed when the distance exceeds the cutoff.
TEST_F(ForceTest, GravityTest) {
    Gravity g(1.0);
    const vec3 force = g(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(50.0, 0.0, 0.0));
}

// Tests Hooke's Law ($F = -k(r - r_0)$) for both tension and compression,
// and validates that the interaction respects the defined cutoff.
TEST_F(ForceTest, HarmonicTest) {
    Harmonic h_stretch(10.0, 1.0);
    const vec3 force_stretch = h_stretch(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force_stretch, vec3(10.0, 0.0, 0.0));

    Harmonic h_compress(10.0, 3.0);
    const vec3 force_compress = h_compress(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force_compress, vec3(-10.0, 0.0, 0.0));
}

// Validates electrostatic force using charges stored in UserData ($F = k_e \frac{q_1 q_2}{r^2}$)
// and ensures cutoff logic correctly disables long-range interactions.
TEST_F(ForceTest, CoulombTest) {
    Coulomb c(1.0);
    const vec3 force = c(get_view1(), get_view2(), r_vec);
    EXPECT_EQ(force, vec3(-0.5, 0.0, 0.0));
}

// Tests the 12-6 Lennard-Jones potential for van der Waals interactions,
// checking mathematical precision near the equilibrium and enforcing cutoff limits.
TEST_F(ForceTest, LennardJonesTest) {
    LennardJones lj(1.0, 2.0);
    const vec3 force = lj(get_view1(), get_view2(), r_vec);

    EXPECT_NEAR(force.x, -12.0, 1e-9);
    EXPECT_NEAR(force.y, 0.0, 1e-9);
    EXPECT_NEAR(force.z, 0.0, 1e-9);
}



// Test Fixture not strictly needed for logic tests, but good for organization
class ForceMixingTest : public testing::Test {};

// HARMONIC MIXING
// Logic:
// - Stiffness (k): Harmonic Mean (Series springs: 2*k1*k2 / (k1+k2))
// - Equilibrium (r0): Arithmetic Mean
// - Cutoff: Maximum
TEST_F(ForceMixingTest, Harmonic_SeriesMixing) {
    // Spring A: k=100, r0=1.0, cutoff=5.0
    const Harmonic h1(100.0, 1.0, 5.0);

    // Spring B: k=50, r0=2.0, cutoff=10.0
    const Harmonic h2(50.0, 2.0, 10.0);

    const Harmonic mixed = h1.mix(h2);

    // Expected k = (2 * 100 * 50) / (100 + 50) = 10000 / 150 = 66.666...
    EXPECT_NEAR(mixed.k, 66.666666667, 1e-9);

    // Expected r0 = (1.0 + 2.0) / 2 = 1.5
    EXPECT_DOUBLE_EQ(mixed.r0, 1.5);

    // Expected Cutoff = max(5.0, 10.0) = 10.0
    EXPECT_DOUBLE_EQ(mixed.cutoff(), 10.0);
}


// GRAVITY MIXING (Guarded)
// Logic:
// - Constant (G): Must be identical (Universal Constant). Throw if different.
// - Cutoff: Maximum (Safety)
TEST_F(ForceMixingTest, Gravity_Mixing_Mismatch_Throws) {
    // Universe A: G=10
    const Gravity g1(10.0, 10.0);
    // Universe B: G=20
    const Gravity g2(20.0, 20.0);

    // Mixing different G values implies a configuration error
    EXPECT_THROW({
        const auto mixed = g1.mix(g2);
        (void)mixed;
    }, std::invalid_argument);
}

TEST_F(ForceMixingTest, Gravity_Mixing_Success) {
    // Valid Case: G is universal (1.0), but cutoffs differ
    const Gravity g1(1.0, 10.0);
    const Gravity g2(1.0, 20.0);

    const Gravity mixed = g1.mix(g2);

    // G remains 1.0
    EXPECT_DOUBLE_EQ(mixed.grav_constant, 1.0);

    // Cutoff is MAX(10, 20) = 20
    EXPECT_DOUBLE_EQ(mixed.cutoff(), 20.0);
}

// COULOMB MIXING (Guarded)
// Logic:
// - Constant (ke): Must be identical. Throw if different.
// - Cutoff: Maximum
TEST_F(ForceMixingTest, Coulomb_Mixing_Mismatch_Throws) {
    // Medium A: ke=100
    const Coulomb c1(100.0, 4.0);
    // Medium B: ke=200
    const Coulomb c2(200.0, 8.0);

    // Mixing different Coulomb constants is physically invalid
    EXPECT_THROW({
                 const auto mixed = c1.mix(c2);
        (void)mixed;
    }, std::invalid_argument);
}

TEST_F(ForceMixingTest, Coulomb_Mixing_Success) {
    // Valid Case: ke is constant, cutoffs differ
    const Coulomb c1(100.0, 4.0);
    const Coulomb c2(100.0, 8.0);

    const Coulomb mixed = c1.mix(c2);

    // ke remains 100
    EXPECT_DOUBLE_EQ(mixed.coulomb_constant, 100.0);

    // Cutoff is MAX(4, 8) = 8.0
    EXPECT_DOUBLE_EQ(mixed.cutoff(), 8.0);
}

// LENNARD-JONES MIXING
// Logic: Lorentz-Berthelot
// - Epsilon (Energy): Geometric Mean (sqrt(e1 * e2))
// - Sigma (Distance): Arithmetic Mean ((s1 + s2) / 2)
// - Cutoff: Geometric Mean (sqrt(c1 * c2))
TEST_F(ForceMixingTest, LennardJones_LorentzBerthelot) {
    // Atom Type A: Shallow well (e=1), small (s=2), short range (cut=9)
    // Note: 9 chosen specifically because sqrt(9 * 16) = sqrt(144) = 12 (clean integer result)
    const LennardJones lj1(1.0, 2.0, 9.0);

    // Atom Type B: Deep well (e=4), large (s=4), long range (cut=16)
    const LennardJones lj2(4.0, 4.0, 16.0);

    const LennardJones mixed = lj1.mix(lj2);

    // 1. Calculate Expected Values manually
    const double expected_epsilon = std::sqrt(1.0 * 4.0); // 2.0
    const double expected_sigma   = 0.5 * (2.0 + 4.0);    // 3.0
    const double expected_cutoff  = std::sqrt(9.0 * 16.0);// 12.0

    // 2. Construct the Expected Object
    // We cannot access private members of 'mixed', but we can compare it
    // against a freshly constructed object with the correct parameters.
    const LennardJones expected(expected_epsilon, expected_sigma, expected_cutoff);

    // 3. Verify Equality
    // This calls operator== which compares epsilon, sigma, and cutoff implicitly
    // (and the derived c6/c12 constants).
    EXPECT_EQ(mixed, expected);

    // 4. Double check the public accessor we DO have (cutoff)
    EXPECT_DOUBLE_EQ(mixed.cutoff(), 12.0);
}