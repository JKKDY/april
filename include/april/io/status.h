#pragma once
#include <iomanip>

#include "april/io/monitor.h"

namespace april::io {
	class ProgressBar : public Monitor {
	public:
		using Monitor::Monitor;

		void record(const size_t step, double, const Particles&) const {
			constexpr size_t bar_width = 50;
			const float progress = static_cast<float>(step + 1) / static_cast<float>(num_steps);
			const auto pos = static_cast<size_t>(bar_width * progress);

			std::cout << "\r[";
			for (size_t i = 0; i < bar_width; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << std::setw(3) << static_cast<int>(progress * 100.0) << "%";
			std::cout.flush();

			if (step + 1 == num_steps) {
				std::cout << std::endl;
			}
		}
	};
}

