#include <april/april.hpp>
#include <filesystem>
#include "april/containers/linked_cells/lc_aos.hpp"
#include "april/containers/linked_cells/lc_soa.hpp"
#include "april/containers/linked_cells/lc_aosoa.hpp"
#include "april/containers/linked_cells/lc_config.hpp"
#include "april/forces/lennard_jones.hpp"


using namespace april;
namespace fs = std::filesystem;



static constexpr int NX = 50, NY = 50, NZ = 50;
static constexpr double a = 1.1225;
static constexpr double sigma = 1.0;
static constexpr double epsilon = 3.0;
static constexpr double r_cut = 3 * sigma;

// Grid physical span
static constexpr double Lx = (NX - 1) * a;
static constexpr double Ly = (NY - 1) * a;
static constexpr double Lz = (NZ - 1) * a;

struct LJNoCutoff : LennardJones {
	using LennardJones::LennardJones;
	// cutoff2 (which should return the cutoff^2 is used for the distance check in the force update function
	// we hijack it to return a large value such that the distance check never succeeds
	// thus causing the system to evaluate the force between every particle pair regardless of actual distance
	[[nodiscard]] double cutoff2() const {
		return 1e100;
	}

	[[nodiscard]] LJNoCutoff mix(LJNoCutoff const&) const noexcept {
		return {0,0,0};
	}
};


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

	for (int i = 0; i < 10; i ++) {

		Environment env (forces<LJNoCutoff>, boundaries<Reflective>);
		env.add_particles(grid);
		env.set_origin(origin);
		env.set_extent(extent);
		env.add_force(LJNoCutoff(epsilon, sigma, r_cut), to_type(0));
		env.set_boundaries(Reflective(), all_faces);

		const auto container = container::LinkedCellsAoS()
			// .with_cell_size(container::CellSize::Cutoff)
			.with_abs_cell_size(3.0)
			.with_cell_ordering(hilbert_order)
			.with_block_size(8);

		auto system = build_system(env, container);


		size_t n_interactions = 0;

		constexpr double dt = 0.000001;
		constexpr int steps  = 25;

		monitor::BenchmarkResult bench;
		VelocityVerlet integrator(system, monitors<Benchmark, ProgressBar>);
		integrator.add_monitor(Benchmark(&bench));
		integrator.run_for_steps(dt, steps);


		n_interactions = 0;
		system.for_each_interaction_pair<+env::Field::position>(
			[&](auto, auto, auto) {
				n_interactions++;
			}
		);
		std::cout << "#interactions: " << n_interactions * steps << std::endl;
		std::cout << "ns / interaction: " << bench.integration_time_s / (n_interactions * steps) * 1e9 << std::endl;
	}
}

