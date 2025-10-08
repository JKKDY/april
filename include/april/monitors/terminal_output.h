#pragma once

#include "april/monitors/monitor.h"
#include "april/env/particle.h"

namespace april::monitor {
	class TerminalOutput final : public Monitor {
	public:
		explicit TerminalOutput(const size_t write_frequency = 1): Monitor(write_frequency) {}

		void record(const size_t step, double, const Particles& particles) {
			std::cerr << "step: " << step <<  "\n";
			for (const auto & p : particles) {
				std::cerr << p.to_string() << "\n";
			}
		}
	};
}