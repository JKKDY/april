
#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

// This file is used for fast testing during developments
// consequently it does not contain a fixed scenario
int main() {

	Force

	// struct ForceTest : force::Force {
	// 	ForceTest() : Force(-1) {}
	//
	// 	[[nodiscard]] vec3 eval(env::internal::Particle const& , env::internal::Particle const& , const vec3 & ) const {
	// 		return {};
	// 	}
	//
	// 	ForceTest mix(const ForceTest &) {
	// 		return ForceTest();
	// 	}
	// };
	//
	// env::internal::Particle p1;
	// env::internal::Particle p2;
	//
	// ForceTest f;
	// f(p1,p2, vec3{});

	// const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/sandbox";
	// remove_all(dir_path);   // delete the directory and all contents
	// create_directory(dir_path); // recreate the empty directory
	//
	// auto drop = ParticleSphere()
	// 	.at({150,150,0})
	// 	.radius_xyz({20, 20, 0})
	// 	.mass(1.0)
	// 	.spacing(1)
	// 	.type(0);
	//
	// // auto thermostat = VelocityScalingThermostat(0.5,
	// // 	controller::TemperatureNotSet,
	// // 	controller::TemperatureNotSet,
	// // 	Trigger::every(1000));
	//
	// auto gravity = UniformField({0, -12.44, 0});
	//
	// auto env = Environment(forces<LennardJones, NoForce>,
	// 	boundaries<Reflective>,
	// 	controllers<VelocityScalingThermostat>,
	// 	fields<UniformField>
	// )
	// 	.with_extent(303,180, 0)
	// 	.with_force(NoForce(), to_type(0))
	// 	.with_particles(drop)
	// 	.with_boundaries(Reflective(), all_faces)
	// 	// .with_controller(thermostat)
	// 	.with_field(gravity);
	//
	// auto container = LinkedCells();
	// auto system = build_system(env, container);
	//
	// auto integrator = StoermerVerlet(system, monitors<Benchmark, ProgressBar, BinaryOutput>)
	// 	.with_monitor(Benchmark())
	// 	.with_monitor(BinaryOutput(Trigger::every(100), dir_path))
	// 	.with_monitor(ProgressBar(Trigger::every(100)))
	// 	.with_dt(0.0002)
	// 	.for_duration(5)
	// 	.run();
}
