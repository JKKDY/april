
#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/two_body_collision";

	ParticleCuboid cuboid1;
	cuboid1.origin = {0,0,0};
	cuboid1.mean_velocity = {0,0,0};
	cuboid1.particle_count = {40, 8, 1};
	cuboid1.mass = 1;
	cuboid1.distance = 1.1225;
	cuboid1.type = 0;

	ParticleCuboid cuboid2;
	cuboid2.origin = {15,15,0};
	cuboid2.mean_velocity = {0,-10,0};
	cuboid2.particle_count = {8, 8, 1};
	cuboid2.mass = 1;
	cuboid2.distance = 1.1225;
	cuboid2.type = 0;

	Environment env;
	env.add_particle_cuboid(cuboid1);
	env.add_particle_cuboid(cuboid2);
	env.set_extent({60,50,1});
	env.set_origin({-10,-10,0});
	env.add_force_to_type(LennardJones(5, 1), 0);

	auto algo = DirectSum();
	auto system = compile(env, algo);

	// using Monitors = april::io::Monitors<BinaryOutput, ProgressBar, Benchmark>;


	auto integrator = StoermerVerlet(system);
	// integrator.add_monitor(BinaryOutput(50, dir_path.string()));
	integrator.add_monitor(ProgressBar(10));
	integrator.add_monitor(Benchmark());
	integrator.run(0.0002, 2);
}
