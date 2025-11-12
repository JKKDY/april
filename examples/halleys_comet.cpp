#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/halleys_comet";
	remove_all(dir_path);   // delete the directory and all contents
	create_directory(dir_path); // recreate the empty directory

	const auto env = Environment (forces<PowerLaw>)
	                 .with_particle({0.0, 0.0, 0.0},     {0.0, 0.0, 0.0},     1.0)
	                 .with_particle({0.0, 1.0, 0.0},     {-1.0, 0.0, 0.0},    3.0e-6)
	                 .with_particle({0.0, 5.36, 0.0},    {-0.425, 0.0, 0.0},  9.55e-4)
	                 .with_particle({34.75, 0.0, 0.0},   {0.0, 0.0296, 0.0},  1.0e-14)
	                 .with_force(PowerLaw(2), to_type(0))
	                 .with_extent({50,50, 0});

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	auto integrator = Yoshida4(system, monitors<BinaryOutput, ProgressBar, Benchmark>)
		.with_monitor(BinaryOutput(Trigger::every(50), dir_path.string()))
		.with_monitor(ProgressBar(Trigger::every(50)))
		.with_monitor(Benchmark())
		.run_for_duration(0.014, 1000);
}

