#include <gtest/gtest.h>
#include <random>

#include "utils.h"
#include "april/april.hpp"
#include "april/containers/linked_cells.hpp"

using namespace april;

// --------------------------------------
// TEST CONFIGURATION & TYPED SUITE SETUP
// --------------------------------------
template <typename ContainerT, typename OrderingT, template<typename...> typename IntegratorT>struct TestConfig {

    static auto create(double cutoff) {
        auto c = ContainerT{};

        // Only apply LinkedCells-specific configuration if the methods exist
        if constexpr (requires { c.with_abs_cell_size(cutoff); }) {
            c.with_abs_cell_size(cutoff)
             .with_block_size(2)
             .with_skin_factor(0.0);
        }

        return OrderingT::apply(std::move(c));
    }

    template<typename S>
    static auto make_integrator(S& sys) {
        return IntegratorT<S>(sys);
    }
};


struct OrderDefault { static auto apply(auto&& c) { return std::forward<decltype(c)>(c); } };
struct OrderHilbert { static auto apply(auto&& c) { return c.with_cell_ordering(hilbert_order); } };

using PhysicsConfigs = testing::Types<
    TestConfig<DirectSum<Layout::AoS>, OrderDefault, VelocityVerlet>,
    TestConfig<DirectSum<Layout::AoSoA<>>, OrderDefault, Yoshida4>,
    TestConfig<LinkedCells<Layout::AoSoA<8>>, OrderHilbert, VelocityVerlet>,
    TestConfig<LinkedCells<Layout::AoSoA<16>>, OrderHilbert, Yoshida4>
>;

template <typename T>
class PhysicsConservationTest : public testing::Test {};
TYPED_TEST_SUITE(PhysicsConservationTest, PhysicsConfigs);



// -----------------
// PHYSICS UTILITIES
// -----------------
namespace physics_test {
    // We must manually calculate potential energy here because april::LennardJones
    // currently only evaluates forces.
    double lj_potential_energy(double r2, double eps, double sig, double rc2) {
        if (r2 > rc2) return 0.0;
        const double s2 = (sig * sig) / r2;
        const double s6 = s2 * s2 * s2;
        return 4.0 * eps * (s6 * s6 - s6);
    }

    template<typename S>
    double calculate_total_energy(const S& system, double eps, double sig, double rc) {
        double kinetic = 0.0;
        double potential = 0.0;
        const double rc2 = rc * rc;

        constexpr auto F = ParticleField::position | ParticleField::velocity | ParticleField::mass | ParticleField::id;

        system.template for_each_particle_view<ParallelPolicy::Serial>(
            april::scalar_kernel<F>([&](auto&& p1) {
                kinetic += 0.5 * p1.mass * p1.velocity.norm_squared();

                // Inner loop for exact O(N^2) potential check
                system.template for_each_particle_view<ParallelPolicy::Serial>(
                    april::scalar_kernel<F>([&](auto&& p2) {
                        // Use IDs to guarantee unique pair evaluation (mimics j = i + 1)
                        if (p1.id < p2.id) {
                            const double r2 = (p2.position - p1.position).norm_squared();
                            potential += lj_potential_energy(r2, eps, sig, rc2);
                        }
                    })
                );
            })
        );
        return kinetic + potential;
    }

    template<typename S>
    vec3 calculate_total_momentum(const S& system) {
        vec3 P = {0, 0, 0};
        constexpr auto F = ParticleField::velocity | ParticleField::mass;

        system.template for_each_particle_view<ParallelPolicy::Serial>(
            april::scalar_kernel<F>([&](auto&& p) {
                P += p.mass * p.velocity;
            })
        );
        return P;
    }
}



// -----
// TESTS
// -----

// checks net force between particles is 0 (newton 3: ∑ F_i = 0)
TYPED_TEST(PhysicsConservationTest, NetForceZero_N3L) {
    constexpr double rc = 3.0;
    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({10, 10, 10});

    // Create jittered grid (4x4x4 = 64 particles)
    for (int i = 0; i < 64; ++i) {
        vec3 pos = {1.0 + (i%4)*2.0, 1.0 + ((i/4)%4)*2.0, 1.0 + (i/16)*2.0};
        env.add_particle(make_particle(0, pos, {0,0,0}, 1.0));
    }
    env.add_force(LennardJones(5.0, 1.0, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    sys.update_forces();

    vec3 net_force = {0, 0, 0};
    sys.template for_each_particle_view<ParallelPolicy::Serial>(
        april::scalar_kernel<ParticleField::force>([&](auto&& p) {
            net_force += p.force;
        })
    );

    // Threshold for 64 particles in double precision accumulation
    EXPECT_NEAR(net_force.x, 0.0, 1e-11);
    EXPECT_NEAR(net_force.y, 0.0, 1e-11);
    EXPECT_NEAR(net_force.z, 0.0, 1e-11);
}

// test that overall momentum is conserved (∑ m_i * v_i =const)
TYPED_TEST(PhysicsConservationTest, MomentumConservation) {
    const double rc = 2.5;
    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({10, 10, 10});

    // Initialize two particles with opposing velocities
    env.add_particle(make_particle(0, {4.5, 5.0, 5.0}, {1.0, 0, 0}, 1.0));
    env.add_particle(make_particle(0, {5.5, 5.0, 5.0}, {-1.0, 0, 0}, 1.0));
    env.add_force(LennardJones(5.0, 1.0, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    const vec3 P0 = physics_test::calculate_total_momentum(sys);

    VelocityVerlet integrator(sys, monitors<>);
    integrator.run_for_steps(0.001, 200);

    const vec3 P_final = physics_test::calculate_total_momentum(sys);

    EXPECT_NEAR(P0.x, P_final.x, 1e-13);
    EXPECT_NEAR(P0.y, P_final.y, 1e-13);
    EXPECT_NEAR(P0.z, P_final.z, 1e-13);
}

// test that energy is conserved
TYPED_TEST(PhysicsConservationTest, HamiltonianEnergyConservation) {
    constexpr double rc = 3.0;
    constexpr double eps = 5.0;
    constexpr double sig = 1.0;
    constexpr double dt = 0.0005;

    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({15, 15, 15});

    // Add 8 particles in a sparse cube to ensure they interact but don't explode
    for (int i = 0; i < 8; ++i) {
        vec3 pos = {5.0 + (i%2)*2.0, 5.0 + ((i/2)%2)*2.0, 5.0 + (i/4)*2.0};
        env.add_particle(make_particle(0, pos, {0.1, -0.1, 0}, 1.0, ParticleState::ALIVE, i));
    }
    env.add_force(LennardJones(eps, sig, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    
    const double E0 = physics_test::calculate_total_energy(sys, eps, sig, rc);
    
    VelocityVerlet integrator(sys, monitors<>);
    integrator.run_for_steps(dt, 500);

    const double E_final = physics_test::calculate_total_energy(sys, eps, sig, rc);
    const double rel_drift = std::abs(E_final - E0) / std::abs(E0);

    // Velocity-Verlet typically stays within 1e-5 for small dt
    EXPECT_LT(rel_drift, 1e-4) << "Energy drift too high! Initial: " << E0 << " Final: " << E_final;
}


// Long-duration stress test to catch rare numerical drifts or race conditions
TYPED_TEST(PhysicsConservationTest, StressTest_EnergyConservation) {
#ifndef NDEBUG
    GTEST_SKIP() << "[ SKIP ] Stress tests require optimization (-O3). Skipping in Debug mode.";
#endif

    // 1. Setup Environment
    constexpr double rc  = 2.5;
    constexpr double eps = 1.0;
    constexpr double sig = 1.0;
    constexpr double dt  = 0.0005;

    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({20.0, 20.0, 20.0});

    // 2. Initialize Dense Fluid (6x6x6 Grid)
    for (int k = 0; k < 6; ++k) {
        for (int j = 0; j < 6; ++j) {
            for (int i = 0; i < 6; ++i) {
                vec3 pos = {2.0 + i * 1.1, 2.0 + j * 1.1, 2.0 + k * 1.1};
                env.add_particle(make_particle(0, pos, {0, 0, 0}, 1.0));
            }
        }
    }
    env.add_force(LennardJones(eps, sig, rc), to_type(0));

    // 3. Build System & Integrator
    auto sys = build_system(env, TypeParam::create(rc));
    VelocityVerlet integrator(sys, monitors<>);

    // 4. Relaxation Phase (Warm-up)
    // Allow the artificial grid to settle before measuring conservation
    integrator.run_for_steps(dt, 100);

    // 5. Measurement Phase
    const double E0 = physics_test::calculate_total_energy(sys, eps, sig, rc);

    integrator.run_for_steps(dt, 5000);

    const double E_final = physics_test::calculate_total_energy(sys, eps, sig, rc);
    const double relative_drift = std::abs(E_final - E0) / std::abs(E0);

    // 6. Assertion
    // Relative drift should stay below 0.5% for 5000 steps in NVE
    EXPECT_LT(relative_drift, 5e-3)
        << "Major energy drift! E_start: " << E0 << " | E_end: " << E_final;
}