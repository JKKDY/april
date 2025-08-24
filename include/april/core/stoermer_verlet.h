#pragma once
#include "april/core/integrator.h"
#include "april/env/particle.h"
#include "april/io/monitor.h"

namespace april::core {
	template <io::IsMonitor ... TMonitors>
	class StoermerVerlet : public impl::Integrator<TMonitors...> {
	public:
		using State = env::ParticleState;
		using Base = impl::Integrator<TMonitors...>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		void integration_step() const {
			for (auto i = sys.index_start(); i <= sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_position(dt * p.velocity + (dt*dt) / (2 * p.mass) * p.force);
			}

			sys.update_forces();

			for (auto i = sys.index_start(); i <= sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_velocity(dt / 2 / p.mass * (p.force + p.old_force));
			}
		}
	};
}