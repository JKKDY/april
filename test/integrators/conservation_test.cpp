#include <gtest/gtest.h>
#include <random>

#include "constant_force.h"
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
        return IntegratorT<S, monitor::internal::MonitorPack<>>(sys);
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
    env.add_interaction(LennardJones(5.0, 1.0, rc), to_type(0));

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
    env.add_interaction(LennardJones(5.0, 1.0, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    const vec3 P0 = physics_test::calculate_total_momentum(sys);

    auto integrator = TypeParam::make_integrator(sys);
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
    env.add_interaction(LennardJones(eps, sig, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    
    const double E0 = physics_test::calculate_total_energy(sys, eps, sig, rc);
    
    auto integrator = TypeParam::make_integrator(sys);
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
    env.add_interaction(LennardJones(eps, sig, rc), to_type(0));

    // 3. Build System & Integrator
    auto sys = build_system(env, TypeParam::create(rc));
    auto integrator = TypeParam::make_integrator(sys);

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


// -------------------------
// ADD THESE TO YOUR TESTS
// -------------------------

// Test that reversing the velocities retraces the path exactly (Symplectic property)
TYPED_TEST(PhysicsConservationTest, TimeReversibility) {
    constexpr double rc = 3.0;
    constexpr double dt = 0.001;
    constexpr int steps = 200;

    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({15, 15, 15});

    // Setup a small, interacting cluster
    env.add_particle(make_particle(0, {7.0, 7.5, 7.5}, {1.0, 0, 0}, 1.0, ParticleState::ALIVE, 0));
    env.add_particle(make_particle(0, {8.0, 7.5, 7.5}, {-1.0, 0, 0}, 1.0, ParticleState::ALIVE, 1));
    env.add_interaction(LennardJones(5.0, 1.0, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));

    // Store initial state
    const vec3 pos0 = get_particle_by_id(sys, 0).position;
    const vec3 vel0 = get_particle_by_id(sys, 0).velocity;

    auto integrator = TypeParam::make_integrator(sys);

    // run Forward
    integrator.run_for_steps(dt, steps);

    // reverse Velocities using scalar_kernel
    sys.template for_each_particle<ParallelPolicy::Serial>(
        april::scalar_kernel<ParticleField::velocity, ParticleField::velocity>(
            [](auto&& p) { p.velocity = -p.velocity; }
        )
    );

    // 4. Run Backward (Forward integration with reversed velocities)
    integrator.run_for_steps(dt, steps);

    // 5. Verification
    const vec3 pos_final = get_particle_by_id(sys, 0).position;
    const vec3 vel_final = get_particle_by_id(sys, 0).velocity;

    // Position should return to start; Velocity should be exactly reversed
    EXPECT_NEAR(pos_final.x, pos0.x, 1e-12);
    EXPECT_NEAR(vel_final.x, -vel0.x, 1e-12);
}

TYPED_TEST(PhysicsConservationTest, RespectsStationaryState) {
    // Using a constant force to ensure significant displacement for ALIVE particles
    Environment env(forces<ConstantForce>);
    env.set_extent({10, 10, 10});

    // P0: Stationary at origin.
    env.add_particle(make_particle(0, {0,0,0}, {0,0,0}, 1.0, ParticleState::STATIONARY, 0));
    // P1: Alive at {5,5,5}.
    env.add_particle(make_particle(0, {5,5,5}, {0,0,0}, 1.0, ParticleState::ALIVE, 1));

    // Push both particles in the +X direction
    env.add_interaction(ConstantForce(10.0, 0.0, 0.0), to_type(0));

    auto sys = build_system(env, TypeParam::create(5.0));
    auto integrator = TypeParam::make_integrator(sys);

    // Run for enough steps to see clear movement
    integrator.run_for_steps(0.01, 20);

    auto p0 = get_particle_by_id(sys, 0);
    auto p1 = get_particle_by_id(sys, 1);

    // 1. Stationary particle must remain exactly at its initial coordinates
    EXPECT_EQ(p0.position, vec3(0,0,0))
        << "Stationary particle moved! Position update mask is failing.";
    EXPECT_EQ(p0.velocity, vec3(0,0,0))
        << "Stationary particle accelerated! Velocity update mask is failing.";

    // 2. Alive particle must have moved significantly
    EXPECT_GT(p1.position.x, 5.01);
    EXPECT_GT(p1.velocity.x, 0.0);
}



// ----------------------------------------------------
// INTEGRATOR CONVERGENCE SETUP (To avoid name clashes)
// ----------------------------------------------------

template <
    typename ContainerT,
    typename OrderingT,
    template<typename...> typename IntegratorT,
    size_t Order // The mathematical order (2 for VV, 4 for Yoshida)
>
struct ConvergenceConfig {
    static constexpr size_t expected_order = Order;

    static auto create(double cutoff) {
        auto c = ContainerT{};
        if constexpr (requires { c.with_abs_cell_size(cutoff); }) {
            c.with_abs_cell_size(cutoff).with_block_size(2).with_skin_factor(0.0);
        }
        return OrderingT::apply(std::move(c));
    }

    template<typename S>
    static auto make_integrator(S& sys) {
        return IntegratorT<S, monitor::internal::MonitorPack<>>(sys);
    }
};

using ConvergenceConfigs = testing::Types<
    ConvergenceConfig<DirectSum<Layout::AoS>, OrderDefault, VelocityVerlet, 2>,
    ConvergenceConfig<DirectSum<Layout::AoSoA<>>, OrderDefault, Yoshida4, 4>,
    ConvergenceConfig<LinkedCells<Layout::AoSoA<8>>, OrderHilbert, VelocityVerlet, 2>,
    ConvergenceConfig<LinkedCells<Layout::AoSoA<16>>, OrderHilbert, Yoshida4, 4>
>;

template <typename T>
class IntegratorConvergenceTest : public testing::Test {};
TYPED_TEST_SUITE(IntegratorConvergenceTest, ConvergenceConfigs);

// Verifies the integrator follows its derived order of accuracy O(dt^n)
TYPED_TEST(IntegratorConvergenceTest, OrderOfAccuracy) {
    constexpr double G = 1.0, M = 1.0, R = 1.0, v = 1.0;
    constexpr double duration = 0.5;

    auto get_orbit_error = [&](double dt) {
        Environment env(forces<Gravity>);
        env.add_particle(make_particle(0, {0,0,0}, {0,0,0}, M, ParticleState::STATIONARY, 0));
        env.add_particle(make_particle(0, {0,R,0}, {v,0,0}, 1e-6, ParticleState::ALIVE, 1));
        env.add_interaction(Gravity(G), to_type(0));
        env.set_extent({10,10,10});

        auto sys = build_system(env, TypeParam::create(R*2));
        auto integrator = TypeParam::make_integrator(sys);

        integrator.run_for_duration(dt, duration);

        const double exact_x = R * std::sin(duration);
        const double exact_y = R * std::cos(duration);
        const vec3 exact_pos(exact_x, exact_y, 0.0);

        // 2. Error is the distance between numerical result and analytical truth
        const vec3 final_pos = get_particle_by_id(sys, 1).position;
        return (final_pos - exact_pos).norm();
    };

    const double dt1 = 0.001;
    const double dt2 = 0.0005;

    const double err1 = get_orbit_error(dt1);
    const double err2 = get_orbit_error(dt2);
    const double ratio = err1 / err2;

    if constexpr (TypeParam::expected_order == 2) {
        // For O(dt^2), halving dt should reduce error by ~4x
        EXPECT_GT(ratio, 3.5) << "2nd order convergence failed. Ratio: " << ratio;
    } else if constexpr (TypeParam::expected_order == 4) {
        // For O(dt^4), halving dt should reduce error by ~16x
        EXPECT_GT(ratio, 12.0) << "4th order convergence failed. Ratio: " << ratio;
    }
}
