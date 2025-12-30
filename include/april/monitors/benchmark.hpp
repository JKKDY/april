#pragma once

#include <chrono>
#include <numeric>
#include <algorithm> // for sort
#include <cmath>     // for sqrt
#include <iomanip>   // for setprecision

#include "april/monitors/monitor.hpp"

namespace april::monitor {
	class Benchmark : public Monitor {
	public:
		Benchmark() : Monitor(shared::Trigger::always()) {}

		void initialize() {
			glob_start_time = std::chrono::high_resolution_clock::now();
			timings.reserve(num_steps);
		}

		template<class S>
		void before_step(const core::SystemContext<S> & sys) {
			updates += sys.size();
			start_time = std::chrono::high_resolution_clock::now();
		}

		template<class S>
		void record(const core::SystemContext<S> &) {
			end_time = std::chrono::high_resolution_clock::now();
			const auto elapsed = std::chrono::duration<double>(end_time - start_time).count();
			timings.push_back(elapsed);
		}

		void finalize() {
			if (timings.empty()) return;

			glob_end_time = std::chrono::high_resolution_clock::now();

			// basic sums
          const double glob_total_s = std::chrono::duration<double>(glob_end_time - glob_start_time).count();
          const double total_integ_s = std::accumulate(timings.begin(), timings.end(), 0.0);
          const size_t steps = timings.size();

          // averages and throughput
          const double avg_s = total_integ_s / static_cast<double>(steps);
          const double its_per_sec = 1.0 / avg_s; // Iterations per second
          const double mups = (static_cast<double>(updates) / total_integ_s) / 1'000'000.0;

          // statistical analysis (Min, Max, Median)
          // we sort a copy to find median/percentiles without altering original order if needed later
          std::vector<double> sorted_timings = timings;
          std::ranges::sort(sorted_timings);

          const double min_s = sorted_timings.front();
          const double max_s = sorted_timings.back();
          const double median_s = sorted_timings[steps / 2];

          // standard deviation
          double variance_accum = 0.0;
          for (const double t : timings) {
             variance_accum += (t - avg_s) * (t - avg_s);
          }
          const double std_dev = std::sqrt(variance_accum / static_cast<double>(steps));

          // Output
          std::cout << "\n" << std::string(40, '-') << "\n";
          std::cout << " [APRIL BENCHMARK REPORT] \n";
          std::cout << std::string(40, '-') << "\n";

          std::cout << std::fixed << std::setprecision(5);
          std::cout << "  Steps processed:    " << steps << "\n";
          std::cout << "  Particles processed: " << updates << "\n";
          std::cout << "  Wall time (total):  " << glob_total_s << " s\n";
          std::cout << "  Integration time:   " << total_integ_s << " s\n";
          std::cout << std::string(40, '-') << "\n";

          std::cout << std::setprecision(2);
          std::cout << "  Throughput:         " << its_per_sec << " it/s\n";
          std::cout << "  Performance:        " << mups << " MUPS\n";
          std::cout << std::string(40, '-') << "\n";

          std::cout << std::setprecision(6);
          std::cout << "  Avg step time:      " << avg_s << " s\n";
          std::cout << "  Median step time:   " << median_s << " s\n";
          std::cout << "  Min step time:      " << min_s << " s\n";
          std::cout << "  Max step time:      " << max_s << " s\n";
          std::cout << "  Std Deviation:      " << std_dev << " s\n";
          std::cout << std::string(40, '-') << "\n\n";
		}

	private:
		using Clock = std::chrono::high_resolution_clock;
		Clock::time_point glob_start_time;
		Clock::time_point glob_end_time;
		Clock::time_point start_time;
		Clock::time_point end_time;
		std::vector<double> timings;
		uint64_t updates = 0;
	};

}

