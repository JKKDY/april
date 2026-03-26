#include <gtest/gtest.h>
#include <random>
#include <numeric>

#include "april/april.hpp"
#include "april/containers/linked_cells.hpp"
#include "april/exec/executors/sequential_executor.hpp"

using namespace april;

// ---------------------------------------------------------
// 1. TEST CONFIGURATION & TYPED SUITE SETUP
// ---------------------------------------------------------

template <typename ContainerT, typename OrderingT>
struct TestConfig {
    static auto create(double cutoff) {
        auto c = ContainerT{}
            .with_abs_cell_size(cutoff)
            .with_block_size(2)
            .with_skin_factor(0.0);
        return OrderingT::apply(std::move(c));
    }
};

struct OrderDefault { static auto apply(auto&& c) { return std::forward<decltype(c)>(c); } };
struct OrderHilbert { static auto apply(auto&& c) { return c.with_cell_ordering(hilbert_order); } };

using PhysicsConfigs = testing::Types<
    TestConfig<LinkedCells<Layout::AoS>, OrderDefault>,
    TestConfig<LinkedCells<Layout::AoSoA<8>>, OrderHilbert>
>;

template <typename T>
class PhysicsConservationTest : public testing::Test {};
TYPED_TEST_SUITE(PhysicsConservationTest, PhysicsConfigs);

// ---------------------------------------------------------
// 2. PHYSICS HELPERS
// ---------------------------------------------------------

namespace physics_test {
    // Standard LJ potential for energy calculation: V(r) = 4 * eps * ((sig/r)^12 - (sig/r)^6)
    double lj_potential(double r2, double eps, double sig, double rc2) {
        if (r2 > rc2) return 0.0;
        double s2 = (sig * sig) / r2;
        double s6 = s2 * s2 * s2;
        double shift = 0.0; // For strict conservation, you might want to shift the potential at rc
        return 4.0 * eps * (s6 * s6 - s6);
    }

    template<typename S>
    double calculate_total_energy(S& system, double eps, double sig, double rc) {
        double kinetic = 0.0;
        double potential = 0.0;
        const double rc2 = rc * rc;
        const size_t n = system.particle_count();

        for (size_t i = 0; i < n; ++i) {
            auto p = system.get_particle(i);
            kinetic += 0.5 * p.mass * p.velocity.norm_squared();
            
            // Potential (O(N^2) here for the 'Golden Truth' check)
            for (size_t j = i + 1; j < n; ++j) {
                auto p2 = system.get_particle(j);
                double r2 = (p2.position - p.position).norm_squared();
                potential += lj_potential(r2, eps, sig, rc2);
            }
        }
        return kinetic + potential;
    }

    template<typename S>
    vec3 calculate_total_momentum(S& system) {
        vec3 P = {0, 0, 0};
        for (size_t i = 0; i < system.particle_count(); ++i) {
            auto p = system.get_particle(i);
            P += p.mass * p.velocity;
        }
        return P;
    }
}

// ---------------------------------------------------------
// 3. THE TESTS
// ---------------------------------------------------------

TYPED_TEST(PhysicsConservationTest, NetForceZero_N3L) {
    const double rc = 3.0;
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
    for (size_t i = 0; i < sys.particle_count(); ++i) {
        net_force += sys.get_particle(i).force;
    }

    // Threshold for 64 particles in double precision
    EXPECT_NEAR(net_force.x, 0.0, 1e-11);
    EXPECT_NEAR(net_force.y, 0.0, 1e-11);
    EXPECT_NEAR(net_force.z, 0.0, 1e-11);
}

TYPED_TEST(PhysicsConservationTest, MomentumConservation) {
    const double rc = 2.5;
    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({10, 10, 10});

    // Initialize two particles with opposing velocities
    env.add_particle(make_particle(0, {4.5, 5, 5}, {1.0, 0, 0}, 1.0));
    env.add_particle(make_particle(0, {5.5, 5, 5}, {-1.0, 0, 0}, 1.0));
    env.add_force(LennardJones(5.0, 1.0, rc), to_type(0));

    auto sys = build_system(env, TypeParam::create(rc));
    const vec3 P0 = physics_test::calculate_total_momentum(sys);

    VelocityVerlet integrator(sys, monitors<>);
    integrator.run_for_steps(0.001, 200);

    const vec3 P_final = physics_test::calculate_test::calculate_total_momentum(sys);

    EXPECT_NEAR(P0.x, P_final.x, 1e-13);
    EXPECT_NEAR(P0.y, P_final.y, 1e-13);
    EXPECT_NEAR(P0.z, P_final.z, 1e-13);
}

TYPED_TEST(PhysicsConservationTest, HamiltonianEnergyConservation) {
    const double rc = 3.0;
    const double eps = 5.0;
    const double sig = 1.0;
    const double dt = 0.0005;

    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.set_extent({15, 15, 15});

    // Add 8 particles in a sparse cube to ensure they interact but don't explode
    for (int i = 0; i < 8; ++i) {
        vec3 pos = {5.0 + (i%2)*2.0, 5.0 + ((i/2)%2)*2.0, 5.0 + (i/4)*2.0};
        env.add_particle(make_particle(0, pos, {0.1, -0.1, 0}, 1.0));
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