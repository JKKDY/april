#pragma once

#include "april/monitors/monitor.h"

namespace april::monitor {
	class TerminalOutput final : public Monitor {
	public:
		using Monitor::Monitor;

		void record(const core::SimulationContext & sys) {
			std::cerr << "step: " << sys.step() <<  "\n";

			for (size_t i = sys.index_start(); i < sys.index_end(); i++) {
				env::ParticleView p = sys.get_particle_by_index(i);
				std::cerr << p.to_string() << "\n";
			}
		}
	};
}