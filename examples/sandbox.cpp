#include <april/april.hpp>
using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	// 1) Generate a block  particles
	auto blob = ParticleCuboid()
		.at(0,0, 10)
		.count(10, 10, 10)
		.spacing(1.2)
		.mass(1.0)
		.type(0)
		.thermal([](vec3 /*pos*/) {
			constexpr double sigma = 1;
			return math::maxwell_boltzmann_velocity(sigma);
		});

	// 2) Define the Environment
	auto env = Environment(
			forces<LennardJones>,
			boundaries<Reflective>,
			fields<UniformField>
		)
		.with_particles(blob)
		.with_extent(30, 30, 50)
		.with_force(LennardJones(3,1), to_type(0))
		.with_field(UniformField({0.0, 0.0, -5}))
		.with_boundaries(Reflective(), all_faces);

	// 3) Build the optimized system (using Linked Cells for O(N) scaling)
	auto container = LinkedCellsAoS();
	auto system = build_system(env, container);

	// 4) Run the simulation
	VelocityVerlet(system, monitors<ProgressBar, BinaryOutput>)
		.with_monitor(ProgressBar(Trigger::every(50)))
		.with_monitor(BinaryOutput(Trigger::every(50), "test/"))
		.with_dt(0.001)
		.for_duration(10)
		.run();
}