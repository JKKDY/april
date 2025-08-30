#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/halleys_comet";
	fs::remove_all(dir_path);   // delete the directory and all contents
	fs::create_directory(dir_path); // recreate the empty directory

	Environment env (forces<InverseSquare>);
	env.add({0.0, 0.0, 0.0},     {0.0, 0.0, 0.0},     1.0);
	env.add({0.0, 1.0, 0.0},     {-1.0, 0.0, 0.0},    3.0e-6);
	env.add({0.0, 5.36, 0.0},    {-0.425, 0.0, 0.0},  9.55e-4);
	env.add({34.75, 0.0, 0.0},   {0.0, 0.0296, 0.0},  1.0e-14);
	env.add_force(InverseSquare(), to_type(0));

	constexpr auto algo = DirectSum();
	auto system = build_system(env, algo);

	StoermerVerlet integrator(system);
	integrator.add_monitor(BinaryOutput(50, dir_path));
	integrator.add_monitor(ProgressBar(10));
	integrator.add_monitor(Benchmark());

	integrator.run_for(0.014, 1000);
}

