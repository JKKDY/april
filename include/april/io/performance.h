#pragma once

#include <chrono>
#include <numeric>

#include "april/io/monitor.h"

namespace april::io {
	class Benchmark : public Monitor {
	public:
		Benchmark() : Monitor(1) {}

		void initialize() {
			glob_start_time = std::chrono::high_resolution_clock::now();
		}

		void before_step(const size_t, double, const Particles&particles) {
			updates+= particles.size();
			start_time = std::chrono::high_resolution_clock::now();
		}

		void record(const size_t, double, const Particles&) {
			end_time = std::chrono::high_resolution_clock::now();
			const auto elapsed = std::chrono::duration<double>(end_time - start_time).count();
			timings.push_back(elapsed);
		}

		void finalize() {
			if (timings.empty()) return;

			glob_end_time = std::chrono::high_resolution_clock::now();
			const double glob_total = std::chrono::duration<double>(glob_end_time - glob_start_time).count();
			const double total = std::accumulate(timings.begin(), timings.end(), 0.0);
			const double avg = total / static_cast<double>(timings.size());
			const double mups = static_cast<double>(updates) / total / 1000000;

			std::cout << "\n[Benchmark Monitor]\n";
			std::cout << "  Total integration time: " << total << " s\n";
			std::cout << "  Total program time:     " << glob_total << " s\n";
			std::cout << "  Avg. per step:          " << avg << " s\n";
			std::cout << "  Avg. MUPS:              " << mups << "MU/s\n";
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

