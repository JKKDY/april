#include <april/april.hpp>
#include <filesystem>
#include <random>
#include <cmath>
#include <vector>

using namespace april;
namespace fs = std::filesystem;

// 1. Define a Custom Force to handle n^2 softening
struct SoftGravity : force::Force {
    double G;
    double eps_sq; // Softening parameter squared

    SoftGravity(const double G, const double eps) : Force(force::no_cutoff), G(G), eps_sq(eps * eps) {}

    // Tell APRIL which particle fields this force needs to access
    static constexpr auto fields = ParticleField::position | ParticleField::mass;

    // The evaluation kernel called by the container
    vec3 eval(auto p1, auto p2, const vec3& r) const noexcept {
        // r is the relative distance vector
        double dist_sq = r.x*r.x + r.y*r.y + r.z*r.z + eps_sq;
        double inv_dist_cube = 1.0 / std::pow(dist_sq, 1.5);
        
        // F = -G * (m1 * m2) / (r^2 + eps^2)^(3/2) * r
        double mag = -G * p1.mass * p2.mass * inv_dist_cube;
        return {r.x * mag, r.y * mag, r.z * mag}; 
    }

    [[nodiscard]] SoftGravity mix(SoftGravity const& other) const {
        if (std::abs(G - other.G) > 1e-9) {
            throw std::invalid_argument("Cannot mix different Gravitational Constants!");
        }
        return SoftGravity(G, std::max(cutoff(), other.cutoff()));
    }
};

int main() {
    const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/spiral_galaxy";
    fs::remove_all(dir_path);   
    fs::create_directory(dir_path); 

    // Simulation Constants
    constexpr int NUM_PARTICLES = 10000;
    constexpr double G_CONST = 1.0;
    constexpr double CORE_MASS = 100000.0;
    constexpr double MAX_RADIUS = 100.0;
    constexpr int NUM_ARMS = 2;
    constexpr double WINDING = 5.0;
    constexpr double THICKNESS = 2.0;

    std::vector<Particle> galaxy;
    galaxy.reserve(NUM_PARTICLES);

    // 2. Generate the Supermassive Core
    galaxy.push_back(Particle()
        .at(0, 0, 0)
        .with_velocity(0, 0, 0)
        .with_mass(CORE_MASS)
        .as_type(0));

    // 3. Generate the Stars
    std::mt19937 gen(42); // Seeded for reproducibility
    std::uniform_real_distribution<double> dist_01(0.0, 1.0);
    std::uniform_real_distribution<double> dist_scatter(-0.5, 0.5);

    for (int i = 0; i < NUM_PARTICLES - 1; ++i) {
        // Density weighted towards the center
        double r = MAX_RADIUS * std::pow(dist_01(gen), 2);
        if (r < 1.0) r = 1.0; // Keep away from exact center

        // Spiral arm math
        int arm = i % NUM_ARMS;
        double offset = (2.0 * M_PI / NUM_ARMS) * arm;
        double theta = offset + (r / MAX_RADIUS) * WINDING + dist_scatter(gen);

        // Position (pinched at the edges)
        double x = r * std::cos(theta);
        double y = r * std::sin(theta);
        double z = dist_scatter(gen) * 2.0 * THICKNESS * (1.0 - r / MAX_RADIUS);

        // Velocity (Stable circular orbit based on the core's mass)
        double v_mag = std::sqrt((G_CONST * CORE_MASS) / r);
        double vx = (-y / r) * v_mag;
        double vy = (x / r) * v_mag;
        double vz = 0.0;

        galaxy.push_back(Particle()
            .at(x, y, z)
            .with_velocity(vx, vy, vz)
            .with_mass(1.0)
            .as_type(0));
    }

    // 4. Configure the APRIL Environment
    // We register our custom SoftGravity force here
    auto env = Environment(forces<SoftGravity>)
        .with_particles(galaxy)
        .with_extent(300, 300, 300)        // Large enough domain to contain the galaxy
        .with_origin(-150, -150, -150)     // Centered at 0,0,0
        // Apply SoftGravity with G=1.0 and epsilon=0.1 to all particles of type 0
        .with_force(SoftGravity(G_CONST, 0.1), to_type(0));

    // 5. Build System using DirectSum (for exact O(n^2) gravity)
    auto container = DirectSum<Layout::AoS>();
    auto system = build_system(env, container);

    // 6. Run the Simulation with Yoshida4
    auto integrator = Yoshida4(system, monitors<ProgressBar, BinaryOutput>)
        .with_monitor(BinaryOutput(Trigger::every(10), dir_path.string()))
        .with_monitor(ProgressBar(Trigger::every(10)))
        .run_for_duration(0.005, 20); // Small timestep, 2000 steps

    return 0;
}