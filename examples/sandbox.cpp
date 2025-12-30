
#include <april/april.hpp>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

// This file is used for fast testing during developments
// consequently it does not contain a fixed scenario
int main() {

	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	const auto env = Environment(forces<Gravity>)
	                 .with_particle(Particle().at(0,0,0).as_type(0).with_mass(1))
	                 .with_particle(Particle().at(1,0,0).as_type(0).with_mass(1))
	                 .with_force(Gravity(), to_type(0));

	auto system = build_system(env, DirectSum());

	auto integrator = VelocityVerlet(system, monitors<BinaryOutput, Benchmark, TerminalOutput>)
		.with_monitor(BinaryOutput(Trigger::always(), dir_path.string()))
		.with_monitor(TerminalOutput(Trigger::always()))
		.run_for_steps(0.001, 1000);
}
