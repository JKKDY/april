#include <iostream>
#include <vector>
#include <chrono>

#include "april/common.hpp"
#include "april/particle/particle.hpp"

// --- Configuration ---
static constexpr int NX = 40, NY = 40, NZ = 40;
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
    vec3 old_position = {};
    vec3 force = {};
    vec3 old_force = {};
    vec3 velocity = {};
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
    std::cout<< sizeof(Particle) << std::endl;
    std::cout<< sizeof(april::env::internal::ParticleRecord<april::env::NoUserData>) << std::endl;

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

    // Simulation Loop
    for (int step = 0; step < STEPS; ++step) {

        for (auto & p : particles) {
            p.old_position = p.position;
            p.position += DT * p.velocity + (DT * DT) / (2 * MASS) * p.force;
        }

        for (auto & p: particles) {
            p.old_force = p.force;
            p.force = {};
        }


        for (size_t i = 0; i < N; ++i) {
            auto & p1 = particles[i];


            for (size_t j = i + 1; j < N; ++j) {
                auto & p2 = particles[j];

                vec3 r = p2.position - p1.position;
                const double r2 = r.norm_squared();
                if (r2 < R_CUT2){
                    const vec3 f = force(r);

                    p1.force += f;
                    p2.force -= f;
                }
            }
        }

        for (auto & p : particles) {
            p.velocity +=  DT / 2 / MASS * (p.force + p.old_force);
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