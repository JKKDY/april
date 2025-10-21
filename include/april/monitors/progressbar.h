#pragma once
#include <iomanip>

#include "april/monitors/monitor.h"

namespace april::monitor {
	class ProgressBar : public Monitor {
	public:
		using Monitor::Monitor;

		void record(const core::SimulationContext & sys) const {
			constexpr size_t bar_width = 50;
			const float progress = static_cast<float>(sys.step() + 1) / static_cast<float>(num_steps);
			const auto pos = static_cast<size_t>(bar_width * progress);

			std::cout << "\r[";
			for (size_t i = 0; i < bar_width; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
			std::cout.flush();

			if (sys.step() + 1 == num_steps) {
				std::cout << std::endl;
			}
		}
	};
}

