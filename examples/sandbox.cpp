
#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

// This file is used for fast testing during developments
// consequently it does not contain a fixed scenario
int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	Environment env (forces<NoForce, LennardJones, InverseSquare>, boundaries<Reflective, Absorb, Outflow, Periodic>);
	env.add_particle({0.0, 0.5, 0.5},     {2, 0.0, 0.0},     1.0);
	env.add_particle({0.5, 0.0, 0.0},     {-1, 0.0, 0.0},     1.0);

	env.add_particle({0.0, 0.0, 0.5},     {0, 1.0, 0.0},     1.0);
	env.add_particle({0.5, 0.5, 0.0},     {0, -2.0, 0.0},     1.0);

	env.add_particle({0.0, 0.0, 0.0},     {0, 0.0, 2.0},     1.0);
	env.add_particle({0.5, 0.5, 0.5},     {0, 0.0, -3.0},     1.0);

	env.set_extent({2,2,2});
	env.set_origin({-1,-1,-1});
	env.add_force(InverseSquare(), to_type(0));
	env.set_boundaries(Periodic(), all_faces);

	auto container = LinkedCells(1);
	auto system = build_system(env, container);

	auto integrator = StoermerVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>);
	integrator.add_monitor(BinaryOutput(50, dir_path));
	integrator.add_monitor(Benchmark());
	integrator.add_monitor(ProgressBar(100));
	integrator.run_for(0.0005, 5);
}
