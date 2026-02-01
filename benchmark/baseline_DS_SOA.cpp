#include <iostream>
#include <vector>
#include <chrono>

// --- Configuration ---
static constexpr int NX = 20;
static constexpr int NY = 20;
static constexpr int NZ = 20;
static constexpr double A = 1.1225;
static constexpr double MASS = 1.0;

// Lennard-Jones Parameters
static constexpr double SIGMA = 1.0;
static constexpr double EPSILON = 5.0;
static constexpr double R_CUT = 3.0 * SIGMA;

// Derived constants
static constexpr double SIGMA2 = SIGMA * SIGMA;
static constexpr double R_CUT2 = R_CUT * R_CUT;

// Simulation settings
static constexpr double DT = 0.0002;
static constexpr int STEPS = 200;
// Pre-calculate constant for half-kick: 0.5 * dt / mass
static constexpr double DT_HALF_MASS = 0.5 * DT / MASS;

int main() {
    // 1. Initialization
    size_t N_reserved = NX * NY * NZ;

    std::vector<double> rx, ry, rz;
    std::vector<double> vx, vy, vz;
    std::vector<double> fx, fy, fz;
    std::vector<double> old_fx, old_fy, old_fz;

    rx.reserve(N_reserved); ry.reserve(N_reserved); rz.reserve(N_reserved);
    vx.reserve(N_reserved); vy.reserve(N_reserved); vz.reserve(N_reserved);
    fx.reserve(N_reserved); fy.reserve(N_reserved); fz.reserve(N_reserved);
    old_fx.reserve(N_reserved); old_fy.reserve(N_reserved); old_fz.reserve(N_reserved);

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
                rx.push_back(i * A + off_x);
                ry.push_back(j * A + off_y);
                rz.push_back(k * A + off_z);

                vx.push_back(0.0); vy.push_back(0.0); vz.push_back(0.0);
                fx.push_back(0.0); fy.push_back(0.0); fz.push_back(0.0);
                old_fx.push_back(0.0); old_fy.push_back(0.0); old_fz.push_back(0.0);
            }
        }
    }

    size_t N = rx.size();

    std::cout << "Starting Benchmark (Kick-Drift-Kick)\n";
    std::cout << "Particles: " << N << "\n";
    std::cout << "Steps: " << STEPS << "\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    // 2. Simulation Loop
    for (int step = 0; step < STEPS; ++step) {

        for (size_t i = 0; i < N; ++i) {
            rx[i] += vx[i] * DT + (DT*DT) * 0.5 / MASS * fx[i];
            ry[i] += vy[i] * DT + (DT*DT) * 0.5 / MASS * fy[i];
            rz[i] += vz[i] * DT + (DT*DT) * 0.5 / MASS * fz[i];

            old_fx[i] = fx[i];
            old_fy[i] = fy[i];
            old_fz[i] = fz[i];

            fx[i] = 0.0;
            fy[i] = 0.0;
            fz[i] = 0.0;
        }

        // force update
        for (size_t i = 0; i < N; ++i) {
            double i_rx = rx[i];
            double i_ry = ry[i];
            double i_rz = rz[i];

            double i_fx = 0.0;
            double i_fy = 0.0;
            double i_fz = 0.0;

            for (size_t j = i + 1; j < N; ++j) {
                double dx = i_rx - rx[j];
                double dy = i_ry - ry[j];
                double dz = i_rz - rz[j];

                double r2 = dx*dx + dy*dy + dz*dz;

                if (r2 < R_CUT2) {
                    double r2inv = 1.0 / r2;
                    double s2 = SIGMA2 * r2inv;
                    double s6 = s2 * s2 * s2;
                    double s12 = s6 * s6;

                    double f_scalar = 24.0 * EPSILON * r2inv * (2.0 * s12 - s6);

                    double f_x = f_scalar * dx;
                    double f_y = f_scalar * dy;
                    double f_z = f_scalar * dz;

                    i_fx += f_x;
                    i_fy += f_y;
                    i_fz += f_z;

                    fx[j] -= f_x;
                    fy[j] -= f_y;
                    fz[j] -= f_z;
                }
            }
            fx[i] += i_fx;
            fy[i] += i_fy;
            fz[i] += i_fz;
        }

        // velocity update
        for (size_t i = 0; i < N; ++i) {
            vx[i] += (fx[i] + old_fx[i]) * DT_HALF_MASS;
            vy[i] += (fy[i] + old_fy[i]) * DT_HALF_MASS;
            vz[i] += (fz[i] + old_fz[i]) * DT_HALF_MASS;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;

    std::cout << "Done.\n";
    std::cout << "Time elapsed: " << diff.count() << " s\n";
    std::cout << "Steps/sec: " << STEPS / diff.count() << "\n";
    std::cout << "Pair interactions/sec: " << (double(N)*double(N)/2.0 * STEPS) / diff.count() << "\n";

    return 0;
}