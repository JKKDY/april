
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



	struct Charge {
		double charge;
	};

	auto cuboid1 = ParticleCuboid{}
	.at({0, 0, 0})
	.velocity({0, 0, 0})
	.count({40, 8, 5})
	.mass(1.0)
	.spacing(1.1225)
	.type(0)
	.with_data(Charge(1));

	auto cuboid2 = ParticleCuboid{}
	.at({15, 15, 0})
	.velocity({0, -20, 0})
	.count({8, 8, 5})
	.mass(1.0)
	.spacing(1.1225)
	.type(0)
	.with_data(Charge(1));


	auto env = Environment(forces<LennardJones>, boundaries<Reflective, Absorb, Open>, particle_data<Charge>)
	   .with_particles(cuboid1)
	   .with_particles(cuboid2)
	   .with_extent(100,80,40)
	   .with_origin(-20,-20,-20)
	   .with_force(LennardJones(5, 1), to_type(0))
	   .with_boundaries(Reflective(), all_faces);

	auto container = LinkedCells();
	auto system = build_system(env, container);

	auto integrator = StoermerVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>)
		.with_monitor(Benchmark())
		.with_monitor(BinaryOutput(Trigger::every(100), dir_path))
		.with_monitor(ProgressBar(Trigger::every(100)))
		.with_dt(0.0002)
		.for_duration(5)
		.run();
}
