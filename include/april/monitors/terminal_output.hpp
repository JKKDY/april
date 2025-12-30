#pragma once

#include <iostream>
#include "april/system/context.hpp"
#include "april/monitors/monitor.hpp"

namespace april::monitor {

	class TerminalOutput final : public Monitor {
	public:
		static constexpr env::FieldMask fields = to_field_mask(env::Field::all);

		using Monitor::Monitor;


		template<class S>
		void record(const core::SystemContext<S> & sys) {
			std::cout << "\n ##########  step: " << sys.step() <<  "  ########## \n";

			for (size_t i = 0; i < sys.size(); ++i) {
				env::ParticleView p = sys.template view<fields>(i);
				std::cout << env::particle_to_string(p) << "\n";
			}
		}
	};
}