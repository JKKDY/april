#include <april/april.h>
#include <filesystem>

using namespace april;
namespace fs = std::filesystem;

int main() {
	Environment env;
	env.add_particle({0.0, 0.0, 0.0},     {0.0, 0.0, 0.0},     1.0);
	env.add_particle({0.0, 1.0, 0.0},     {-1.0, 0.0, 0.0},    3.0e-6);
	env.add_particle({0.0, 5.36, 0.0},    {-0.425, 0.0, 0.0},  9.55e-4);
	env.add_particle({34.75, 0.0, 0.0},   {0.0, 0.0296, 0.0},  1.0e-14);
	env.add_force_to_type(InverseSquare(), 0);

	BinaryOutput out(50, fs::path(PROJECT_SOURCE_DIR) / "output/halleys_comet");

	StoermerVerlet<BinaryOutput> integrator(env);
	integrator.add_monitor(out);
	integrator.run(0.014, 1000);
}

