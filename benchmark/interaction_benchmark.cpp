#include <april/april.hpp>
#include <filesystem>
#include "april/containers/linked_cells/lc_aos.hpp"
#include "april/containers/linked_cells/lc_soa.hpp"
#include "april/containers/linked_cells/lc_aosoa.hpp"
#include "april/containers/linked_cells/lc_config.hpp"

using namespace april;
namespace fs = std::filesystem;



// class InteractionCounter : public Monitor {
// public:
// 	using Monitor::Monitor;
//
// 	size_t n_interactions;
//
// 	template<class S>
// 	void record(const SystemContext<S> & sys) {
//
// 	}
//
// 	void finalize() const {
// 		std::cout << 2 * n_interactions << std::endl;
// 	}
// };



static constexpr int NX = 100, NY = 100, NZ = 100;
static constexpr double a = 1.1225;
static constexpr double sigma = 1.0;
static constexpr double epsilon = 3.0;
static constexpr double r_cut = 3*3 * sigma;

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

	const auto container = container::LinkedCellsSoA()
		// .with_cell_size(container::CellSize::Cutoff)
		.with_abs_cell_size(3.0)
		.with_cell_ordering(hilbert_order)
		.with_block_size(8);

	auto system = build_system(env, container);


	size_t n_interactions = 0;
	size_t pairs_looked_at = 0;

	constexpr double dt = 0.000001;
	constexpr int steps  = 1;

	VelocityVerlet integrator(system, monitors<Benchmark, ProgressBar>);
	integrator.add_monitor(Benchmark());
	integrator.run_for_steps(dt, steps);


	n_interactions = 0;
	system.for_each_interaction_pair<+env::Field::position>(
		[&](auto, auto, const vec3 & r) {
			if (r.norm() < r_cut) n_interactions++;
			pairs_looked_at++;
		}
	);
	std::cout << "#interactions" <<n_interactions * steps << std::endl;
	std::cout << "#pairs looked at during iteration:" << pairs_looked_at * steps << std::endl;
	std::cout << "#particles updated: " + std::to_string(n_interactions*steps*2) << std::endl; // *2 because of newton 3


}

