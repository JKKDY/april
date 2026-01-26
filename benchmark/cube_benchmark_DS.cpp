#include <april/april.hpp>
#include <filesystem>


using namespace april;
namespace fs = std::filesystem;


static constexpr int NX = 20, NY = 20, NZ = 20;
static constexpr double a = 1.1225;
static constexpr double sigma = 1.0;
static constexpr double epsilon = 3.0;
static constexpr double r_cut = 3 * sigma;

// Grid physical span
static constexpr double Lx = (NX - 1) * a;
static constexpr double Ly = (NY - 1) * a;
static constexpr double Lz = (NZ - 1) * a;


int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/bench";
	const vec3 box = {Lx, Ly, Lz};

	ParticleCuboid grid = ParticleCuboid{}
		.at(-0.5*box)
		.velocity({0,0,0})
		.count({NX, NY, NZ})
		.mass(1.0)
		.spacing(a)
		.type(0);

	// Box with margin >= r_cut around grid (non-periodic)
	const vec3 extent = 1.5 * box;
	const vec3 origin = - 0.5 * extent;

	Environment env (forces<LennardJones>, boundaries<Reflective>);
	env.add_particles(grid);
	env.set_origin(origin);
	env.set_extent(extent);
	env.add_force(LennardJones(epsilon, sigma, r_cut), to_type(0));
	env.set_boundaries(Reflective(), all_faces);

	constexpr auto container = DirectSumSoA();
	auto system = build_system(env, container);

	constexpr double dt = 0.0002;
	constexpr int steps  = 100;

	VelocityVerlet integrator(system, monitors<Benchmark>);
	integrator.add_monitor(Benchmark());
	// integrator.add_monitor(ProgressBar(Trigger::every(200)));
	integrator.run_for_steps(dt, steps);

	std::cout << "Particles: " << NX * NY * NZ << "\n"
			 << "Steps: " << steps << "\n"
			 << "dt: " << dt << "\n";
}

