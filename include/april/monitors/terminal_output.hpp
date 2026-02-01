#pragma once

#include <string>
#include <sstream>

#include <iostream>
#include "april/system/context.hpp"
#include "april/monitors/monitor.hpp"

namespace april::monitor {

	class TerminalOutput final : public Monitor {
	public:
		static constexpr env::FieldMask fields = +env::Field::all;

		using Monitor::Monitor;

		template<typename P>
		static std::string particle_to_string(const P & p) {
			std::ostringstream oss;
			oss << "Particle ID: " << p.id << "\n"
				<< "Position: " << p.position.to_string() << "\n"
				<< "Velocity: " << p.velocity.to_string() << "\n"
				<< "Force: " << p.force.to_string() << "\n"
				<< "Mass: " << p.mass << "\n"
				<< "Type: " << p.type << "\n"
				<< "State: " << static_cast<int>(p.state) << "\n";
			return oss.str();
		}


		template<class S>
		void record(const core::SystemContext<S> & sys) {
			std::cout << "\n ##########  step: " << sys.step() <<  "  ########## \n";

			sys.template for_each_particle_view<fields>([&](const auto & p) {
				std::cout << particle_to_string(p) << "\n";
			});
		}
	};
}