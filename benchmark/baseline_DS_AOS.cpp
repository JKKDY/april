#include <iostream>
#include <vector>
#include <chrono>

#include "april/common.h"

// --- Configuration ---
static constexpr int NX = 40;
static constexpr int NY = 40;
static constexpr int NZ = 40;
static constexpr double A = 1.1225;
double MASS = 1.0;

// Lennard-Jones Parameters
double SIGMA = 1.0;
double EPSILON = 5.0;
double R_CUT = 3 * SIGMA;

// Derived constants
double SIGMA2 = SIGMA * SIGMA;
double R_CUT2 = R_CUT * R_CUT;

// Simulation settings
double DT = 0.0002;
int STEPS = 5;


using april::vec3;

struct Particle {
    vec3 position = {};
    vec3 force = {};
    vec3 velocity = {};
    vec3 old_position = {};
};


 vec3 force(const vec3 & r) {
    const double r2 = r.norm_squared();

    const double r2inv = 1.0 / r2;
    const double s2 = SIGMA2 * r2inv;
    const double s6 = s2 * s2 * s2;
    const double s12 = s6 * s6;

    const double magnitude = 24.0 * EPSILON * r2inv * (2.0 * s12 - s6);

    return - magnitude * r;
}


int main() {
    // std::cout<< sizeof(Particle) << std::endl;
    // std::cout<< sizeof(april::env::internal::ParticleRecord<april::env::NoUserData>) << std::endl;

    // 1. Initialization
    size_t N = NX * NY * NZ;

    std::vector<Particle> particles;
    particles.reserve(N);

    // Grid Setup
    double Lx = (NX - 1) * A;
    double Ly = (NY - 1) * A;
    double Lz = (NZ - 1) * A;

    double off_x = -0.5 * Lx;
    double off_y = -0.5 * Ly;
    double off_z = -0.5 * Lz;

    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
            for (int i = 0; i < NX; ++i) {
                auto p = Particle();
                p.position = {i * A + off_x, j * A + off_y, k * A + off_z};
                particles.push_back(p);
            }
        }
    }


    std::cout << "Starting Benchmark (Kick-Drift-Kick)\n";
    std::cout << "Particles: " << N << "\n";
    std::cout << "Steps: " << STEPS << "\n";

    const auto start_time = std::chrono::high_resolution_clock::now();

     const int TILE = 1024;
     Particle* __restrict parts = particles.data();

    // Simulation Loop
    for (int step = 0; step < STEPS; ++step) {
        // 1. KICK + DRIFT (Standard)
        for (auto & p : particles) {
            p.old_position = p.position; // Actually unused in force loop, but part of struct
            p.velocity += (DT * 0.5 / MASS) * p.force;
            p.position += DT * p.velocity;
            p.force = {0.0, 0.0, 0.0}; // Reset forces
        }

        // 2. FORCE (Optimized)
        // Tiled Force Calculation
        for (size_t i_block = 0; i_block < N; i_block += TILE) {
            // End index for outer block
            size_t i_end = std::min((size_t)N, i_block + TILE);

            for (size_t j_block = i_block; j_block < N; j_block += TILE) {
                // End index for inner block
                size_t j_end = std::min((size_t)N, j_block + TILE);

                // Now we process the interactions between these two small blocks
                for (size_t i = i_block; i < i_end; ++i) {

                    // Optimization from before: Load p1 ONCE
                    vec3 p1_pos = parts[i].position;
                    vec3 p1_force_acc = {0.0, 0.0, 0.0};

                    // Careful with the start index to avoid double counting
                    // If blocks are same, start at i+1. If different, start at beginning of j_block.
                    size_t j_start = (i_block == j_block) ? i + 1 : j_block;

                    for (size_t j = j_start; j < j_end; ++j) {

                        vec3 r = parts[j].position - p1_pos;

                        // Early Rejection (Manhattan) - Cheap filter
                        if (std::abs(r.x) > R_CUT || std::abs(r.y) > R_CUT || std::abs(r.z) > R_CUT)
                            continue;

                        double r2 = r.norm_squared();
                        if (r2 < R_CUT2) {
                            // OPTIMIZATION 2: Manual Inline + Common Subexpression Elimination
                            double r2inv = 1.0 / r2;
                            double s2 = SIGMA2 * r2inv;
                            double s6 = s2 * s2 * s2;
                            double s12 = s6 * s6;

                            // Factor out terms.
                            // Original: 24 * eps * r2inv * (2*s12 - s6)
                            double magnitude = 24.0 * EPSILON * r2inv * (2.0 * s12 - s6);

                            vec3 f = -magnitude * r;

                            // OPTIMIZATION 3: Accumulate p1 force in REGISTER
                            p1_force_acc += f;

                            // We still have to write to p2 memory (random access write),
                            // but we've eliminated the write for p1.
                            parts[j].force -= f;
                        }
                    }
                    parts[i].force += p1_force_acc;
                }
            }
        }

        // 3. KICK (Standard)
        for (auto & p : particles) {
            p.velocity += (DT * 0.5 / MASS) * p.force;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
//    std::cout << end_time;
    std::chrono::duration<double> diff = end_time - start_time;

    std::cout << "Done.\n";
    std::cout << "Time elapsed: " << diff.count() << " s\n";
    std::cout << "Steps/sec: " << STEPS / diff.count() << "\n";
    std::cout << "Pair interactions/sec: " << (double(N)*double(N)/2.0 * STEPS) / diff.count() << "\n";

    return 0;
}