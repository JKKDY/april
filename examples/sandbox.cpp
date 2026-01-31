
#include <april/april.hpp>
#include "april/containers/linked_cells/lc_aosoa.hpp"
// #include <filesystem>

using namespace april;
// namespace fs = std::filesystem;

// This file is used for fast testing during developments
// consequently it does not contain a fixed scenario
int main() {

	// const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
	// remove_all(dir_path);   // delete the directory and all contents
	// create_directory(dir_path); // recreate the empty directory
	//
	auto env = Environment(forces<Gravity>)
		// .with_particle(Particle().at(0,0,0).as_type(0).with_mass(1))
		// .with_particle(Particle().at(1,0,0).as_type(1).with_mass(1))
		.with_force(Gravity(), to_type(0))
		.with_force(Gravity(), to_type(1))
		.with_auto_domain(1);

	for (int i = 0; i < 9; i++) {
		env.add_particle(Particle().at(0,0,0).as_type(0).with_mass(1));
		env.add_particle(Particle().at(0,0,0).as_type(1).with_mass(1));
	}
	//
	auto system = build_system(env, container::LinkedCellsAoSoA());
	size_t i = system.capacity();
	std::cout << i << std::endl;
	//
	// auto integrator = VelocityVerlet(system, monitors<BinaryOutput, Benchmark, TerminalOutput>)
	// 	.with_monitor(BinaryOutput(Trigger::always(), dir_path.string()))
	// 	.with_monitor(TerminalOutput(Trigger::always()))
	// 	.run_for_steps(0.001, 1000);



}
