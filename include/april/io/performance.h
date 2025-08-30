#pragma once

#include <chrono>
#include <numeric>

#include "april/io/monitor.h"

namespace april::io {
	class Benchmark : public Monitor {
	public:
		Benchmark() : Monitor(1) {}

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

			const double total = std::accumulate(timings.begin(), timings.end(), 0.0);
			const double avg = total / static_cast<double>(timings.size());
			const double mups = static_cast<double>(updates) / total / 1000000;

			std::cout << "\n[Benchmark Monitor]\n";
			std::cout << "  Total time:   " << total << " s\n";
			std::cout << "  Avg. per step: " << avg << " s\n";
			std::cout << "  Avg. MUPS (mega updates / s): " << mups << "MU/s\n";
		}

	private:
		using Clock = std::chrono::high_resolution_clock;
		Clock::time_point start_time;
		Clock::time_point end_time;
		std::vector<double> timings;
		uint64_t updates = 0;
	};

}

