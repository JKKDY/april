#include <iostream>
#include <vector>
#include <immintrin.h>
#include <chrono>
#include <algorithm>
#include <numeric>



// --- Configuration ---
static constexpr int NX = 40;
static constexpr int NY = 40;
static constexpr int NZ = 40;
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
static constexpr int STEPS = 5;
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

        const double POS_COEFF_VAL = (DT * DT) * 0.5 / MASS;
        const __m256d v_dt       = _mm256_set1_pd(DT);
        const __m256d v_pos_coeff = _mm256_set1_pd(POS_COEFF_VAL);
        const __m256d v_zero     = _mm256_setzero_pd();

        size_t i = 0;
        // Process 4 particles at a time
        for (; i < (N & ~3); i += 4) {
            // --- X COMPONENT ---
            __m256d rx_vec = _mm256_loadu_pd(&rx[i]);
            __m256d vx_vec = _mm256_loadu_pd(&vx[i]);
            __m256d fx_vec = _mm256_loadu_pd(&fx[i]);

            // rx += vx * DT
            rx_vec = _mm256_fmadd_pd(vx_vec, v_dt, rx_vec);
            // rx += fx * COEFF
            rx_vec = _mm256_fmadd_pd(fx_vec, v_pos_coeff, rx_vec);

            _mm256_storeu_pd(&rx[i], rx_vec);     // Update Position
            _mm256_storeu_pd(&old_fx[i], fx_vec); // Copy Force to Old Force
            _mm256_storeu_pd(&fx[i], v_zero);     // Reset Force to 0.0

            // --- Y COMPONENT ---
            __m256d ry_vec = _mm256_loadu_pd(&ry[i]);
            __m256d vy_vec = _mm256_loadu_pd(&vy[i]);
            __m256d fy_vec = _mm256_loadu_pd(&fy[i]);

            ry_vec = _mm256_fmadd_pd(vy_vec, v_dt, ry_vec);
            ry_vec = _mm256_fmadd_pd(fy_vec, v_pos_coeff, ry_vec);

            _mm256_storeu_pd(&ry[i], ry_vec);
            _mm256_storeu_pd(&old_fy[i], fy_vec);
            _mm256_storeu_pd(&fy[i], v_zero);

            // --- Z COMPONENT ---
            __m256d rz_vec = _mm256_loadu_pd(&rz[i]);
            __m256d vz_vec = _mm256_loadu_pd(&vz[i]);
            __m256d fz_vec = _mm256_loadu_pd(&fz[i]);

            rz_vec = _mm256_fmadd_pd(vz_vec, v_dt, rz_vec);
            rz_vec = _mm256_fmadd_pd(fz_vec, v_pos_coeff, rz_vec);

            _mm256_storeu_pd(&rz[i], rz_vec);
            _mm256_storeu_pd(&old_fz[i], fz_vec);
            _mm256_storeu_pd(&fz[i], v_zero);
        }

        // Cleanup Tail (Scalar fallback)
        for (; i < N; ++i) {
            rx[i] += vx[i] * DT + POS_COEFF_VAL * fx[i];
            ry[i] += vy[i] * DT + POS_COEFF_VAL * fy[i];
            rz[i] += vz[i] * DT + POS_COEFF_VAL * fz[i];

            old_fx[i] = fx[i]; old_fy[i] = fy[i]; old_fz[i] = fz[i];
            fx[i] = 0.0;       fy[i] = 0.0;       fz[i] = 0.0;
        }

       // --- FORCE UPDATE (Optimized AVX2 + Early Exit) ---

        // 1. Reset Forces


        // 2. AVX Constants
        const __m256d v_rcut2   = _mm256_set1_pd(R_CUT2);
        const __m256d v_sigma2  = _mm256_set1_pd(SIGMA2);
        const __m256d v_eps24   = _mm256_set1_pd(24.0 * EPSILON);
        const __m256d v_one     = _mm256_set1_pd(1.0);
        const __m256d v_two     = _mm256_set1_pd(2.0);

        // --- OUTER LOOP ---
        for (size_t i = 0; i < N; ++i) {
            __m256d v_ix = _mm256_set1_pd(rx[i]);
            __m256d v_iy = _mm256_set1_pd(ry[i]);
            __m256d v_iz = _mm256_set1_pd(rz[i]);

            // TWO sets of accumulators to break dependency chains
            __m256d v_fx1 = _mm256_setzero_pd(); __m256d v_fx2 = _mm256_setzero_pd();
            __m256d v_fy1 = _mm256_setzero_pd(); __m256d v_fy2 = _mm256_setzero_pd();
            __m256d v_fz1 = _mm256_setzero_pd(); __m256d v_fz2 = _mm256_setzero_pd();

            size_t j = 0;

            // Process 8 particles per iteration (2 chunks of 4)
            for (; j < (N & ~7); j += 8) {

                // --- CHUNK 1 (j to j+3) ---
                __m256d v_jx1 = _mm256_loadu_pd(&rx[j]);
                __m256d v_jy1 = _mm256_loadu_pd(&ry[j]);
                __m256d v_jz1 = _mm256_loadu_pd(&rz[j]);

                __m256d dx1 = _mm256_sub_pd(v_ix, v_jx1);
                __m256d dy1 = _mm256_sub_pd(v_iy, v_jy1);
                __m256d dz1 = _mm256_sub_pd(v_iz, v_jz1);

                __m256d r2_1 = _mm256_fmadd_pd(dx1, dx1, _mm256_fmadd_pd(dy1, dy1, _mm256_mul_pd(dz1, dz1)));

                // --- CHUNK 2 (j+4 to j+7) ---
                __m256d v_jx2 = _mm256_loadu_pd(&rx[j+4]);
                __m256d v_jy2 = _mm256_loadu_pd(&ry[j+4]);
                __m256d v_jz2 = _mm256_loadu_pd(&rz[j+4]);

                __m256d dx2 = _mm256_sub_pd(v_ix, v_jx2);
                __m256d dy2 = _mm256_sub_pd(v_iy, v_jy2);
                __m256d dz2 = _mm256_sub_pd(v_iz, v_jz2);

                __m256d r2_2 = _mm256_fmadd_pd(dx2, dx2, _mm256_fmadd_pd(dy2, dy2, _mm256_mul_pd(dz2, dz2)));

                // --- MASKS ---
                __m256d mask1 = _mm256_cmp_pd(r2_1, v_rcut2, _CMP_LT_OQ);
                __m256d mask2 = _mm256_cmp_pd(r2_2, v_rcut2, _CMP_LT_OQ);

                // Early Exit: Only skip if BOTH chunks are empty
                int m1 = _mm256_movemask_pd(mask1);
                int m2 = _mm256_movemask_pd(mask2);
                if ((m1 | m2) == 0) continue;

                // --- MATH 1 ---
                if (m1) {
                     mask1 = _mm256_and_pd(mask1, _mm256_cmp_pd(r2_1, v_zero, _CMP_GT_OQ));
                     __m256d r2_safe = _mm256_blendv_pd(v_one, r2_1, mask1);
                     __m256d r2inv = _mm256_div_pd(v_one, r2_safe);
                     __m256d s2 = _mm256_mul_pd(v_sigma2, r2inv);
                     __m256d s6 = _mm256_mul_pd(s2, _mm256_mul_pd(s2, s2));
                     __m256d term = _mm256_fmsub_pd(v_two, _mm256_mul_pd(s6, s6), s6);
                     __m256d scalar = _mm256_and_pd(_mm256_mul_pd(_mm256_mul_pd(v_eps24, r2inv), term), mask1);

                     // Accumulate into Set 1
                     v_fx1 = _mm256_fmadd_pd(scalar, dx1, v_fx1);
                     v_fy1 = _mm256_fmadd_pd(scalar, dy1, v_fy1);
                     v_fz1 = _mm256_fmadd_pd(scalar, dz1, v_fz1);
                }

                // --- MATH 2 (Interleaved execution) ---
                if (m2) {
                     mask2 = _mm256_and_pd(mask2, _mm256_cmp_pd(r2_2, v_zero, _CMP_GT_OQ));
                     __m256d r2_safe = _mm256_blendv_pd(v_one, r2_2, mask2);
                     __m256d r2inv = _mm256_div_pd(v_one, r2_safe);
                     __m256d s2 = _mm256_mul_pd(v_sigma2, r2inv);
                     __m256d s6 = _mm256_mul_pd(s2, _mm256_mul_pd(s2, s2));
                     __m256d term = _mm256_fmsub_pd(v_two, _mm256_mul_pd(s6, s6), s6);
                     __m256d scalar = _mm256_and_pd(_mm256_mul_pd(_mm256_mul_pd(v_eps24, r2inv), term), mask2);

                     // Accumulate into Set 2
                     v_fx2 = _mm256_fmadd_pd(scalar, dx2, v_fx2);
                     v_fy2 = _mm256_fmadd_pd(scalar, dy2, v_fy2);
                     v_fz2 = _mm256_fmadd_pd(scalar, dz2, v_fz2);
                }
            }

            // --- MERGE & REDUCE ---
            // Combine the two accumulators
            v_fx1 = _mm256_add_pd(v_fx1, v_fx2);
            v_fy1 = _mm256_add_pd(v_fy1, v_fy2);
            v_fz1 = _mm256_add_pd(v_fz1, v_fz2);

            // Horizontal Sum (same as before)
            double buf_x[4], buf_y[4], buf_z[4];
            _mm256_storeu_pd(buf_x, v_fx1);
            _mm256_storeu_pd(buf_y, v_fy1);
            _mm256_storeu_pd(buf_z, v_fz1);

            fx[i] += buf_x[0] + buf_x[1] + buf_x[2] + buf_x[3];
            fy[i] += buf_y[0] + buf_y[1] + buf_y[2] + buf_y[3];
            fz[i] += buf_z[0] + buf_z[1] + buf_z[2] + buf_z[3];

                    // 6. Tail Loop (Handle N % 4 != 0)
                    for (; j < N; ++j) {
                        double dx = rx[i] - rx[j];
                        double dy = ry[i] - ry[j];
                        double dz = rz[i] - rz[j];
                        double r2 = dx*dx + dy*dy + dz*dz;

                        if (r2 < R_CUT2 && r2 > 0.0) { // Explicitly check > 0 for self-interaction
                            double r2inv = 1.0 / r2;
                            double s2 = SIGMA2 * r2inv;
                            double s6 = s2 * s2 * s2;
                            double s12 = s6 * s6;
                            double f_scalar = 24.0 * EPSILON * r2inv * (2.0 * s12 - s6);

                            fx[i] += f_scalar * dx;
                            fy[i] += f_scalar * dy;
                            fz[i] += f_scalar * dz;
                        }
                    }
                }

        // velocity update
        const __m256d v_dt_half_mass = _mm256_set1_pd(DT_HALF_MASS);

        size_t k = 0; // using 'k' to avoid conflict if you copy-paste in same scope
        for (; k < (N & ~3); k += 4) {
            // --- X ---
            __m256d vx_v = _mm256_loadu_pd(&vx[k]);
            __m256d fx_v = _mm256_loadu_pd(&fx[k]);
            __m256d old_fx_v = _mm256_loadu_pd(&old_fx[k]);

            // sum_forces = fx + old_fx
            __m256d sum_fx = _mm256_add_pd(fx_v, old_fx_v);
            // vx += sum_forces * constant
            vx_v = _mm256_fmadd_pd(sum_fx, v_dt_half_mass, vx_v);
            _mm256_storeu_pd(&vx[k], vx_v);

            // --- Y ---
            __m256d vy_v = _mm256_loadu_pd(&vy[k]);
            __m256d fy_v = _mm256_loadu_pd(&fy[k]);
            __m256d old_fy_v = _mm256_loadu_pd(&old_fy[k]);

            __m256d sum_fy = _mm256_add_pd(fy_v, old_fy_v);
            vy_v = _mm256_fmadd_pd(sum_fy, v_dt_half_mass, vy_v);
            _mm256_storeu_pd(&vy[k], vy_v);

            // --- Z ---
            __m256d vz_v = _mm256_loadu_pd(&vz[k]);
            __m256d fz_v = _mm256_loadu_pd(&fz[k]);
            __m256d old_fz_v = _mm256_loadu_pd(&old_fz[k]);

            __m256d sum_fz = _mm256_add_pd(fz_v, old_fz_v);
            vz_v = _mm256_fmadd_pd(sum_fz, v_dt_half_mass, vz_v);
            _mm256_storeu_pd(&vz[k], vz_v);
        }

        // Cleanup Tail
        for (; k < N; ++k) {
            vx[k] += (fx[k] + old_fx[k]) * DT_HALF_MASS;
            vy[k] += (fy[k] + old_fy[k]) * DT_HALF_MASS;
            vz[k] += (fz[k] + old_fz[k]) * DT_HALF_MASS;
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