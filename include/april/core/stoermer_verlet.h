#pragma once
#include "april/core/integrator.h"
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/io/monitor.h"
#include "april/defaults.h"

namespace april::core {

	template<IsSystem Sys, class Pack> class StoermerVerlet;

	template <IsSystem Sys, class ... TMonitors>
	class StoermerVerlet<Sys, io::MonitorPack<TMonitors...>>
		: public impl::Integrator<Sys, io::MonitorPack<TMonitors...>> {
	public:
		using State = env::ParticleState;
		using Base = impl::Integrator<Sys, io::MonitorPack<TMonitors...>>;
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

	template<class Sys, class... Ms>
	StoermerVerlet(Sys&, io::MonitorPack<Ms...>)
		-> StoermerVerlet<Sys, io::MonitorPack<Ms...>>;

	using core::StoermerVerlet;
	template<class Sys> StoermerVerlet(Sys&)
		-> StoermerVerlet<Sys, DefaultMonitors>;
}