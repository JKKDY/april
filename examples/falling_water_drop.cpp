#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/falling_water_drop";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	auto liquid = ParticleCuboid()
		.at({1.5, 2, 0})
		.velocity({0, 0, 0})
		.count({250, 50, 1})
		.mass(1.0)
		.spacing(1.2)
		.type(0);

	auto drop = ParticleSphere()
		.at({150,150,0})
		.radius_xyz({20, 20, 0})
		.mass(1.0)
		.spacing(1)
		.type(0);

	auto thermostat = VelocityScalingThermostat(0.5,
		controller::TemperatureNotSet,
		controller::TemperatureNotSet,
		Trigger::every(1000));

	auto gravity = UniformField({0, -12.44, 0});

	auto env = Environment(forces<LennardJones>,
		boundaries<Reflective>,
		controllers<VelocityScalingThermostat>,
		fields<UniformField>
	)
		.with_extent(303,180, 0)
		.with_force(LennardJones(1, 1.2), to_type(0))
		.with_particles(liquid)
		.with_particles(drop)
		.with_boundaries(Reflective(), all_faces)
		.with_controller(thermostat)
		.with_field(gravity);

	auto container = LinkedCells();
	auto system = build_system(env, container);

	auto integrator = StoermerVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>)
		.with_monitor(Benchmark())
		.with_monitor(BinaryOutput(Trigger::every(100), dir_path))
		.with_monitor(ProgressBar(Trigger::every(100)))
		.run_for(0.0002, 20);
}
