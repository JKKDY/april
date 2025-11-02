#pragma once

#include <iostream>
#include "april/core/context.h"
#include "april/monitors/monitor.h"

namespace april::monitor {

	template<trigger::IsTrigger Trig>
	class TerminalOutput final : public Monitor<Trig> {
	public:
		static constexpr env::FieldMask fields = to_field_mask(env::Field::all);

		using Monitor<Trig>::Monitor;

		template<class S>
		void record(const core::SystemContext<S> & sys) {
			std::cout << "step: " << sys.step() <<  "\n";

			for (size_t i = sys.index_start(); i < sys.index_end(); ++i) {
				env::ParticleView p = sys.template get_particle_by_index<fields>(i);
				std::cout << p.to_string() << "\n";
			}
		}
	};
}