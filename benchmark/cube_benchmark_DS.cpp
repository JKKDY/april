#include <april/april.h>
#include <filesystem>


using namespace april;
namespace fs = std::filesystem;

// use to check if any particles leave the simulation domain
class ExitMonitor final : public ext::Monitor {
public:
	explicit ExitMonitor(const vec3 & origin_, const vec3 & extent_):
		Monitor(100), extent(extent_), origin(origin_)
	{}

	void record(size_t, double, const Particles& particles) const {
		for (const auto & p : particles) {
			if (
				p.position.x < origin.x or p.position.x > origin.x + extent.x ||
				p.position.y < origin.y or p.position.y > origin.y + extent.y ||
				p.position.z < origin.z or p.position.z > origin.z + extent.z
			) {
				std::cout << "Particle " + std::to_string(p.id) + " at " + p.position.to_string() + " has left domain!";
			}
		}
	}

private:
	vec3 extent;
	vec3 origin;
};


static constexpr int NX = 10, NY = 10, NZ = 10;
static constexpr double a = 1.1225;
static constexpr double sigma = 1.0;
static constexpr double epsilon = 5.0;

// Grid physical span
static constexpr double Lx = (NX - 1) * a;
static constexpr double Ly = (NY - 1) * a;
static constexpr double Lz = (NZ - 1) * a;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/bench";
	const vec3 box = {Lx, Ly, Lz};

	ParticleCuboid grid = ParticleCuboid{}
		.at(-0.5*box)
		.velocity({0,0,0})
		.count({NX, NY, NZ})
		.mass(1.0)
		.spacing(a)
		.type(0);

	// Box with margin >= r_cut around grid (non-periodic)
	const vec3 extent = 1.5 * box;
	const vec3 origin = - 0.5 * extent;

	for (int i = 0; i < 5; i++) {
		Environment env (force<LennardJones>);
		env.add_particles(grid);
		env.set_origin(origin);
		env.set_extent(extent);
		env.add_force(LennardJones(epsilon, sigma, -1), to_type(0));

		constexpr auto container = DirectSum();
		auto system = build_system(env, container);

		constexpr double dt = 0.0002;
		constexpr int steps  = 10000;

		StoermerVerlet integrator(system, monitors<Benchmark>);
		integrator.add_monitor(Benchmark());

		integrator.run_steps(dt, steps);
	}
}

