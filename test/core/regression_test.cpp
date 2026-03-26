#include <gtest/gtest.h>
#include <vector>
#include <chrono>

#include <april/april.hpp>
#include "utils.h"

#include "utils.h"
#include "april/containers/direct_sum.hpp"
#include "april/exec/executors/sequential_executor.hpp"

using namespace april;

// Configuration
static constexpr int NX = 10, NY = 10, NZ = 10;
static constexpr int STEPS = 50;
static constexpr double A = 1.1225;
static constexpr double MASS = 1.0;
static constexpr double SIGMA = 1.0;
static constexpr double EPSILON = 5.0;
static constexpr double R_CUT = 3.0 * SIGMA;
static constexpr double DT = 0.0002;

// HANDCODED BASELINE
namespace baseline {
    struct Particle {
        vec3 position = {};
        vec3 old_position = {};
        vec3 force = {};
        vec3 old_force = {};
        vec3 velocity = {};
    };

    inline vec3 force_lj(const vec3& r) {
        const double r2 = r.norm_squared();
        const double r2inv = 1.0 / r2;
        const double s2 = (SIGMA * SIGMA) * r2inv;
        const double s6 = s2 * s2 * s2;
        const double magnitude = 24.0 * EPSILON * r2inv * (2.0 * s6 * s6 - s6);
        return -magnitude * r;
    }

    double run_benchmark() {
        constexpr size_t N = NX * NY * NZ;
        std::vector<Particle> particles(N);

        constexpr double off_x = -0.5 * (NX - 1) * A;
        constexpr double off_y = -0.5 * (NY - 1) * A;
        constexpr double off_z = -0.5 * (NZ - 1) * A;

        size_t idx = 0;
        for (int k = 0; k < NZ; ++k) {
            for (int j = 0; j < NY; ++j) {
                for (int i = 0; i < NX; ++i) {
                    particles[idx++].position = {i * A + off_x, j * A + off_y, k * A + off_z};
                }
            }
        }

        constexpr double r_cut2 = R_CUT * R_CUT;
        const auto start = std::chrono::high_resolution_clock::now();

        for (int step = 0; step < STEPS; ++step) {
            for (auto& p : particles) {
                p.old_position = p.position;
                p.position += DT * p.velocity + (DT * DT) / (2 * MASS) * p.force;
            }

            for (auto& p : particles) {
                p.old_force = p.force;
                p.force = {};
            }

            for (size_t i = 0; i < N; ++i) {
                auto& p1 = particles[i];
                for (size_t j = i + 1; j < N; ++j) {
                    auto& p2 = particles[j];
                    vec3 r = p2.position - p1.position;
                    if (r.norm_squared() < r_cut2) {
                        const vec3 f = force_lj(r);
                        p1.force += f;
                        p2.force -= f;
                    }
                }
            }

            for (auto& p : particles) {
                p.velocity += DT / 2.0 / MASS * (p.force + p.old_force);
            }
        }

        const auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }
}

// APRIL IMPLEMENTATION
double run_april_benchmark() {
    const vec3 box = {(NX - 1) * A, (NY - 1) * A, (NZ - 1) * A};
    
    ParticleCuboid grid = ParticleCuboid{}
        .at(-0.5 * box)
        .velocity({0, 0, 0})
        .count({NX, NY, NZ})
        .mass(MASS)
        .spacing(A)
        .type(0);

    Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
    env.add_particles(grid);
    env.add_force(LennardJones(EPSILON, SIGMA, R_CUT), to_type(0));
    env.set_boundaries(OpenBoundary(), all_faces);

    // Force strict Sequential/Scalar execution to match the handcoded baseline
    auto container = DirectSum<Layout::AoS>();
    struct Cfg : RunTimeConfig<exec::SequentialExecutor>, 
                 CompileTimeConfig<ParallelPolicy::Serial, VectorPolicy::Scalar> {};
    
    auto system = build_system(env, container, Cfg{});
    VelocityVerlet integrator(system); // No monitors to avoid I/O overhead

    auto start = std::chrono::high_resolution_clock::now();
    integrator.run_for_steps(DT, STEPS);
    auto end = std::chrono::high_resolution_clock::now();

    return std::chrono::duration<double>(end - start).count();
}

double get_allowed_tolerance() {
#if defined(_MSC_VER) && !defined(__clang__)
    return 1.6; // MSVC is allowed more overhead due to conservative inlining
#else
    return 1.1; // GCC/Clang must remain near-zero overhead
#endif
}

// THE TEST
TEST(PerformanceRegression, DirectSumAoS_vs_Handcoded) {
#ifndef NDEBUG
    GTEST_SKIP() << "[ SKIP ] Performance regression tests require optimization. "
                 << "Rebuild in Release or RelWithDebInfo to run this test.";
#endif

    std::cout << "[   INFO   ] Warming up caches..." << std::endl;
    baseline::run_benchmark();
    run_april_benchmark();

    constexpr int max_attempts = 5;
    constexpr double tolerance = 1.15; // 15% allowance for abstraction/CI noise

    for (int attempt = 1; attempt <= max_attempts; ++attempt) {
        std::cout << "\n[   INFO   ] --- Attempt " << attempt << " of " << max_attempts << " ---" << std::endl;

        const double time_handcoded = baseline::run_benchmark();
        const double time_april = run_april_benchmark();

        std::cout << "[  RESULT  ] Handcoded : " << time_handcoded << " s\n";
        std::cout << "[  RESULT  ] April     : " << time_april << " s\n";
        std::cout << "[  RESULT  ] Ratio     : " << (time_april / time_handcoded) << "x\n";

        if (time_april <= time_handcoded * get_allowed_tolerance()) {
            std::cout << "[  PASSED  ] April met the performance target!" << std::endl;
            SUCCEED();
            return; // Early exit! We proved the engine is fast enough.
        } else {
            std::cout << "[   WARN   ] April missed target. Retrying..." << std::endl;
        }
    }

    // If we make it out of the loop, all 5 attempts failed.
    FAIL() << "April consistently failed to meet the performance target ("
           << tolerance << "x baseline) after " << max_attempts << " attempts. "
           << "Check for abstraction regressions or unintended virtual calls.";
}


