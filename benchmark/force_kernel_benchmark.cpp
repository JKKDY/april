#include <april/april.hpp>
#include <filesystem>

#include "april/forces/force.hpp"
#include "april/containers/direct_sum/ds_aos.hpp"
#include "april/containers/direct_sum/ds_soa.hpp"
#include "april/containers/direct_sum/ds_aosoa.hpp"

using namespace april;
namespace fs = std::filesystem;


//###################################################
// IMPORTANT: TO GET COMPARABLE RESULTS BETWEEN APRIL
// AND HARDCODED LOOPS; COMPILE WITH -fno-tree-vectorize
// OTHERWISE SOME OF THE EXAMPLES WILL BE AUTO
// VECTORIZED WHILE APRIL WILL NOT
//###################################################


struct LJ{
	LJ(const double epsilon_, const double sigma_)
	: epsilon(epsilon_), sigma(sigma_) {
		const vec3::type sigma2 = sigma * sigma;
		const vec3::type sigma6 = sigma2 * sigma2 * sigma2;
		const vec3::type sigma12 = sigma6 * sigma6;

		c6_force = 24.0 * epsilon * sigma6;
		c12_force = 48.0 * epsilon * sigma12;
	}

	[[nodiscard]] vec3 eval(const vec3& r) const noexcept {

		const vec3::type inv_r2 = static_cast<vec3::type>(1.0) / (r.x*r.x + r.y*r.y + r.z * r.z);
		const vec3::type inv_r6 = inv_r2 * inv_r2 * inv_r2;
		const vec3::type magnitude = (c12_force * inv_r6 - c6_force) * inv_r6 * inv_r2;

		// Force vector pointing along -r
		return -magnitude * r;
	}
private:
	// Precomputed force constants
	vec3::type c12_force;
	vec3::type c6_force;

	double epsilon; // Depth of the potential well
	double sigma; // Distance at which potential is zero
};


struct LJ_AoS{
	LJ_AoS(const double epsilon_, const double sigma_)
	: epsilon(epsilon_), sigma(sigma_) {
		const vec3::type sigma2 = sigma * sigma;
		const vec3::type sigma6 = sigma2 * sigma2 * sigma2;
		const vec3::type sigma12 = sigma6 * sigma6;

		c6_force = 24.0 * epsilon * sigma6;
		c12_force = 48.0 * epsilon * sigma12;
	}

	[[nodiscard]] auto eval(double x, double y, double z) const noexcept {

		const vec3::type inv_r2 = static_cast<vec3::type>(1.0) / (x*x + y*y + z * z);
		const vec3::type inv_r6 = inv_r2 * inv_r2 * inv_r2;
		const vec3::type magnitude = (c12_force * inv_r6 - c6_force) * inv_r6 * inv_r2;

		// Force vector pointing along -r
		return std::array{-magnitude * x,  -magnitude * y,  -magnitude * z};
	}
private:
	// Precomputed force constants
	vec3::type c12_force;
	vec3::type c6_force;

	double epsilon; // Depth of the potential well
	double sigma; // Distance at which potential is zero
};


constexpr size_t N = 4000;
constexpr double sigma = 1.0;
constexpr double epsilon = 3.0;
constexpr double dt = 0.00001;
constexpr int steps  = 200;
constexpr size_t n_interactions = N * (N+1) / 2 * steps;

int main() {
	const auto dir_path = fs::path(PROJECT_SOURCE_DIR) / "output/bench";

	auto force = LennardJones(epsilon, sigma, force::no_cutoff);

	{
		Environment env (forces<LennardJones>, boundaries<Reflective>);
		env.add_force(force, to_type(0));
		env.set_boundaries(Reflective(), all_faces);

		for (size_t i = 0; i < N; i++) {
			auto p = Particle().at(static_cast<double>(i),0,0).with_mass(1);
			env.add_particle(p);
		}

		constexpr auto container = DirectSumAoS();
		auto system = build_system(env, container);

		monitor::BenchmarkResult bench_results{};
		VelocityVerlet integrator(system, monitors<Benchmark>);
		integrator.add_monitor(Benchmark(&bench_results));
		integrator.run_for_steps(dt, steps);

		const double seconds_per_interaction = bench_results.integration_time_s / static_cast<double>(n_interactions);
		const double ns_per_interaction = seconds_per_interaction * 1e9;

		// std::cout << n_interactions << std::endl;
		// std::cout << bench_results.integration_time_s << std::endl;
		std::cout << "ns/interaction April (full)" << ns_per_interaction << "\n" << std::endl;
	}

	{
		Environment env (forces<LennardJones>, boundaries<Reflective>);
		env.add_force(force, to_type(0));
		env.set_boundaries(Reflective(), all_faces);

		for (size_t i = 0; i < N; i++) {
			auto p = Particle().at(static_cast<double>(i),0,0).with_mass(1);
			env.add_particle(p);
		}

		constexpr auto container = DirectSumAoS();
		auto system = build_system(env, container);

		double total_f_time = 0.0;

		auto start_f = std::chrono::steady_clock::now();
		for (size_t i = 1; i < steps+1; i++) {
			system.update_forces();
		}
		auto end_f = std::chrono::steady_clock::now();

		total_f_time += std::chrono::duration<double>(end_f - start_f).count();
		// std::cout << total_f_time << "\n"<<std::endl;
		std::cout << total_f_time / n_interactions * 1e9 << std::endl;
		std::cout << "ns/interaction April (update_forces call only) " << total_f_time / n_interactions * 1e9 << "\n" << std::endl;
	}



	{
		auto lj = LJ(epsilon, sigma);

		vec3 acc = {};
		double total_f_time = 0.0;

		auto start_f = std::chrono::steady_clock::now();
		for (size_t i = 1; i < n_interactions+1; i++) {
			const auto j = static_cast<double>(i);
			vec3 f = lj.eval(vec3{j} );

#if defined(__GNUC__) || defined(__clang__)
			asm volatile("" : : "r,m"(f) : "memory");
#else
			 acc+= f;
#endif
		}
		auto end_f = std::chrono::steady_clock::now();

		total_f_time += std::chrono::duration<double>(end_f - start_f).count();
		std::cout << "ns/interaction absolute max perf (no memory calls): " << total_f_time / n_interactions * 1e9 << std::endl;

		std::cout << acc.to_string() << std::endl; // printing acc so compiler does not optimize it awway

	}




	{
		auto lj = LJ(epsilon, sigma);

		vec3 acc = {};
		double total_f_time = 0.0;

		std::vector<vec3> pos;
		size_t n = 10000000;

		for (size_t i = 1; i < n+1; i++) {
			pos.emplace_back(0.0001*static_cast<double>(i));
		}

		auto start_f = std::chrono::steady_clock::now();
		for (size_t i = 0; i < n; i++) {
			vec3 f = lj.eval(pos[i]);
			acc+= f;
		}
		auto end_f = std::chrono::steady_clock::now();

		total_f_time += std::chrono::duration<double>(end_f - start_f).count();
		std::cout << "ns/interaction realistic1 (reads from a vector): " << total_f_time / static_cast<double>(n) * 1e9 << std::endl;
		// std::cout << total_f_time << std::endl;

		std::cout << acc.to_string() << std::endl;

	}




	{
		auto lj = LJ(epsilon, sigma);
		double total_f_time = 0.0;

		std::vector<vec3> pos1;
		// std::vector<vec3> pos2;

		size_t n = N;
		std::vector<vec3> f1(n);

		for (size_t i = 1; i < n+1; i++) {
			pos1.emplace_back(0.0001*static_cast<double>(i));
		}

		auto start_f = std::chrono::steady_clock::now();
		for (int k = 0; k < steps; k++) {
			for (size_t i = 0; i < n; i++) {
				vec3 a {};
				for (size_t j = i+1; j < n; j++) {
					vec3 f = lj.eval(pos1[i] - pos1[j]);
					a+= f;
					f1[j]-= f;
				}
				f1[i]+=a;
			}
		}
		auto end_f = std::chrono::steady_clock::now();

		total_f_time += std::chrono::duration<double>(end_f - start_f).count();
		std::cout << "ns/interaction realistic2 (performs triangle traversal & force updates): " << total_f_time / n_interactions * 1e9 << std::endl;
		// std::cout << acc.to_string();
	}


	{
		auto lj = LJ_AoS(epsilon, sigma);

		double total_f_time = 0.0;

		std::vector<double> posx;
		std::vector<double> posy;
		std::vector<double> posz;


		size_t n = N;
		std::vector<double> fx(n);
		std::vector<double> fy(n);
		std::vector<double> fz(n);


		for (size_t i = 1; i < n+1; i++) {
			posx.emplace_back(0.0001*static_cast<double>(i));
			posy.emplace_back(0.0001*static_cast<double>(i));
			posz.emplace_back(0.0001*static_cast<double>(i));

		}

		auto start_f = std::chrono::steady_clock::now();
		for (int k = 0; k < steps; k++) {
			for (size_t i = 0; i < n; i++) {
				double ax = 0,ay = 0,az = 0;
				for (size_t j = i+1; j < n; j++) {
					auto [x,y,z] = lj.eval(posx[i]- posx[j], posy[i]-posy[j], posz[i]-posz[j]);
					ax+= x;
					ay+= y;
					az+= z;

					fx[j]-= x;
					fy[j]-= y;
					fz[j]-= z;

				}
				fx[i]-= ax;
				fy[i]-= ay;
				fz[i]-= az;
			}
		}
		auto end_f = std::chrono::steady_clock::now();

		total_f_time += std::chrono::duration<double>(end_f - start_f).count();
		std::cout << "ns/interaction realistic2 with SoA: " << total_f_time / n_interactions * 1e9 << std::endl;
		// std::cout << acc.to_string();
	}


}
