#pragma once

#include <string>
#include <sstream>

#include <iostream>
#include "april/core/context.hpp"
#include "april/monitors/monitor.hpp"
#include "april/exec/particle_kernel.hpp"

namespace april {

	class TerminalOutput final : public monitor::Monitor {
	public:
		static constexpr auto fields = ParticleField::all;

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

			sys.for_each_particle_view(scalar_kernel<fields>(
				[&](const auto & p) {
				std::cout << particle_to_string(p) << "\n";
			}));
		}
	};
}














