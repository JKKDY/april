
#include <april/april.hpp>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/two_body_collision";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	auto cuboid1 = ParticleCuboid{}
		.at({0, 0, 0})
		.velocity({0, 0, 0})
		.count({40, 8, 5})
		.mass(1.0)
		.spacing(1.1225)
		.type(0);

	auto cuboid2 = ParticleCuboid{}
		.at({15, 15, 0})
		.velocity({0, -20, 0})
		.count({8, 8, 5})
		.mass(1.0)
		.spacing(1.1225)
		.type(0);

	auto env = Environment(forces<LennardJones>, boundaries<Reflective>)
	   .with_particles(cuboid1)
	   .with_particles(cuboid2)
	   .with_extent(100,80,40)
	   .with_origin(-20,-20,-20)
	   .with_force(LennardJones(5, 1), to_type(0))
	   .with_boundaries(Reflective(), all_faces);

	auto container = LinkedCellsAoS();
	auto system = build_system(env, container);

	auto integrator = VelocityVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>)
		.with_monitor(Benchmark())
		.with_monitor(BinaryOutput(Trigger::every(100), dir_path.string()))
		.with_monitor(ProgressBar(Trigger::every(100)))
		.run_for_duration(0.0002, 5);
}



