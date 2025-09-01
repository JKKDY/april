
#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/two_body_collision";
	fs::remove_all(dir_path);   // delete the directory and all contents
	fs::create_directory(dir_path); // recreate the empty directory

	auto cuboid1 = ParticleCuboid{}
		.at({0, 0, 0})
		.velocity({0, 0, 0})
		.count({40, 8, 3})
		.mass(1.0)
		.spacing(1.1225)
		.type(0);

	auto cuboid2 = ParticleCuboid{}
		.at({15, 15, 0})
		.velocity({0, -10, 0})
		.count({8, 8, 3})
		.mass(1.0)
		.spacing(1.1225)
		.type(0);

	Environment env (forces<LennardJones, Harmonic, InverseSquare>);
	env.add(cuboid1);
	env.add(cuboid2);
	env.set_extent({60,50,10});
	env.set_origin({-10,-10,-5});
	env.add_force(LennardJones(5, 1), to_type(0));

	auto container = LinkedCells2();
	auto system = build_system(env, container);

	auto integrator = StoermerVerlet(system, io::monitors<Benchmark, ProgressBar, BinaryOutput>);
	// integrator.add_monitor(BinaryOutput(50, dir_path));
	integrator.add_monitor(Benchmark());
	integrator.add_monitor(ProgressBar(100));
	integrator.run_for(0.0002, 5);
}
