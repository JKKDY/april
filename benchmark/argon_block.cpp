#include <april/april.hpp>
#include <cmath>
#include <iostream>

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace april;

/**
 * @brief Agnostic simulation runner.
 * @tparam Env Environment
 * @tparam Container Container Config
 * @tparam ExecConfig Execution Config
 */
template<
    typename Env,
    typename Container,
    typename ExecConfig>

Benchmark::BenchmarkResult run_simulation(
    Env&& env,
    Container&& container,
    ExecConfig && exec_config,
    size_t warmup_steps,
    size_t bench_steps,
    double dt = 0.001,
    const bool enable_output = false,
    const fs::path& output_path = ""
) {
    auto system = build_system(std::forward<Env>(env), std::forward<Container>(container), exec_config);

    // Warm-up (Silent)
    VelocityVerlet warmup_integrator(system, monitors<>);
    warmup_integrator.run_for_steps(dt, warmup_steps);

    // Benchmark Preparation
    Benchmark::BenchmarkResult res;
    // We include all potential monitors in the static pack
    VelocityVerlet bench_integrator(system, monitors<Benchmark, BinaryOutput, ProgressBar>);

    // Always add the benchmark collector
    bench_integrator.add_monitor(Benchmark(&res));

    // Runtime toggle for expensive IO/Visualization
    if (enable_output && !output_path.empty()) {
        if (!fs::exists(output_path)) fs::create_directories(output_path);
        fs::remove_all(output_path);
        bench_integrator.add_monitor(BinaryOutput(Trigger::every(50), output_path.string()));
        bench_integrator.add_monitor(ProgressBar(Trigger::every(bench_steps / 100)));
    }

    // 5. Execution
    bench_integrator.run_for_steps(dt, bench_steps);

    return res;
}


/**
 * @brief Appends benchmark metadata and results to a CSV file.
 */
inline void save_bench_to_csv(const fs::path& csv_path, const std::string& label, size_t threads, const Benchmark::BenchmarkResult& res) {
    const bool is_new = !fs::exists(csv_path);
    std::ofstream csv(csv_path, std::ios::app);

    if (is_new) {
        csv << "Label,Particles,Threads,Steps,Integration_s,MUPS,Avg_Step_s,Median_Step_s,StdDev_s\n";
    }

    csv << label << ","
        << res.total_updates / res.steps << ","
        << threads << ","
        << res.steps << ","
        << res.integration_time_s << ","
        << res.mups << ","
        << res.avg_step_sec << ","
        << res.median_step_sec << ","
        << res.std_dev_sec << "\n";
}


inline fs::path get_next_run_directory(const fs::path& base_dir) {
    if (!fs::exists(base_dir)) {
        fs::create_directories(base_dir);
        return base_dir / "1";
    }

    int max_id = 0;
    for (const auto& entry : fs::directory_iterator(base_dir)) {
        if (entry.is_directory()) {
            try {
                // Attempt to convert folder name to integer
                int current_id = std::stoi(entry.path().filename().string());
                max_id = std::max(max_id, current_id);
            } catch (...) {}
        }
    }

    return base_dir / std::to_string(max_id + 1);
}

using namespace april;

/**
 * @brief Generates a bulk liquid argon benchmark environment.
 * * @param n_dim Number of particles along one dimension (Total = n_dim^3).
 * @param temperature Initial reduced temperature (T*). Default 1.0.
 * @return A fully configured, declarative Environment.
 */
auto create_argon_environment(const size_t n_dim, const double temperature = 1.0) {
    // Standard LJ liquid parameters
    constexpr double density = 0.8442;
    constexpr double epsilon = 1.0;
    constexpr double sigma = 1.0;
    constexpr double r_cut = 2.5 * sigma;
    constexpr int TYPE_ARGON = 0;

    // Calculate Box Dimensions
    const size_t total_particles = n_dim * n_dim * n_dim;
    const double volume = static_cast<double>(total_particles) / density;
    double L = std::cbrt(volume);
    const double spacing = L / static_cast<double>(n_dim);
    L = L + spacing/2;

    // Center the domain at (0,0,0)
    const vec3 origin = {-L / 2.0, -L / 2.0, -L / 2.0};
    const auto extent = vec3{L, L, L};

    // Velocity generator for thermalization
    auto thermalize = [temperature](vec3 /*pos*/) {
        return math::maxwell_boltzmann_velocity(temperature);
    };

    // Generate the particle grid
    // We offset the start position by spacing/2 to avoid placing particles
    // exactly on the periodic boundary edges.
    auto particle_grid = ParticleCuboid()
        .at(origin + vec3(spacing / 2.0))
        .velocity(vec3(0))
        .count({n_dim, n_dim, n_dim})
        .mass(1.0)
        .spacing(spacing)
        .type(TYPE_ARGON)
        .thermal(thermalize);

    // Build the declarative environment
    auto env = Environment(forces<LennardJones>, boundaries<PeriodicBoundary>)
        .with_extent(extent)
        .with_particles(particle_grid)
        .with_force(LennardJones(epsilon, sigma, r_cut), to_type(TYPE_ARGON))
        .with_boundaries(PeriodicBoundary(), all_faces);

    return env;
}


/**
 * @brief Core benchmarking suite for the Argon Block.
 * Sweeps: Density x System Size x Threads
 */
void run_argon_bench_suite() {
    namespace fs = std::filesystem;
    const auto base_dir = std::filesystem::path(PROJECT_SOURCE_DIR) / "bench_results/argon_sweep";
    const auto master_csv = base_dir / "master_results.csv";
    std::filesystem::create_directories(base_dir);

    // --- Parameter Ranges ---
    const std::vector<double> densities = {0.4, 0.8442, 1.1}; // Gas, Liquid, Near-Solid
    const std::vector<size_t> sizes     = {40, 60, 80, 100, 140, 170, 200}; // 64k to 8M particles
    const std::vector<size_t> threads   = {1, 2, 4, 8, 16, 32, 64}; // Scaling range

    std::cout << ">>> Starting Argon Sweep Suite\n";

    for (double rho : densities) {
        for (size_t n : sizes) {
            for (size_t t : threads) {

                // 1. Prepare unique directory and label
                const fs::path run_path = get_next_run_directory(base_dir);
                fs::create_directories(run_path);

                std::string label = "rho" + std::to_string(rho).substr(0, 4) +
                                    "_n" + std::to_string(n) +
                                    "_t" + std::to_string(t);

                std::cout << "\n[RUN " << run_path.filename().string() << "] " << label << std::endl;

                // 2. Setup Config
                ExecutionConfig cfg;
                cfg.executer_config.n_threads = t;

                auto env = create_argon_environment(n, rho);
                auto container = LinkedCells<Layout::AoSoA<>>().with_absolute_skin(0.3);

                // 3. Execute
                auto result = run_simulation(
                    std::move(env),
                    std::move(container),
                    cfg,
                    25, // warmup
                    500, // bench
                    0.001,
                    false, // No binary output during benchmark
                    run_path
                );

                // 4. Persistence
                save_bench_to_csv(master_csv, label, t, result);
                save_bench_to_csv(run_path / "result.csv", label, t, result);
            }
        }
    }
}


int main() {
    // run_argon_bench_suite();

    ExecutionConfig cfg;
    cfg.executer_config.n_threads = 6;

    auto env = create_argon_environment(100, 1);
    auto container = LinkedCells<Layout::AoSoA<>>().with_absolute_skin(0.3);
    // auto container = DirectSum();

    run_simulation(
      std::move(env),
      std::move(container),
      cfg,
      25, // warmup
      500, // bench
      0.001,
      true, // No binary output during benchmark
      "/mnt/d/Dev/april/animation/output"
  );
}