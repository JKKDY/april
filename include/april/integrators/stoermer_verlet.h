#pragma once
#include "april/integrators/integrator.h"
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/monitors/monitor.h"
#include "april/defaults.h"

namespace april::integrator {

	template<core::IsSystem Sys, class Pack> class StoermerVerlet;

	template <core::IsSystem Sys, class ... TMonitors>
	class StoermerVerlet<Sys, monitor::MonitorPack<TMonitors...>>
		: public impl::Integrator<Sys, monitor::MonitorPack<TMonitors...>> {
	public:
		using State = env::ParticleState;
		using Base = impl::Integrator<Sys, monitor::MonitorPack<TMonitors...>>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		void integration_step() const {
			for (auto i = sys.index_start(); i < sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_position(dt * p.velocity + (dt*dt) / (2 * p.mass) * p.force);
			}

			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			for (auto i = sys.index_start(); i < sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_velocity(dt / 2 / p.mass * (p.force + p.old_force));
			}

			sys.apply_controllers();
		}
	};

	// Deduction guide so user can write StoermerVerlet(sys, MonitorPack<M1, M2, M3>)
	template<class Sys, class... Ms>
	StoermerVerlet(Sys&, monitor::MonitorPack<Ms...>)
		-> StoermerVerlet<Sys, monitor::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys)
	template<class Sys>
	StoermerVerlet(Sys&)
		-> StoermerVerlet<Sys, DefaultMonitors>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<class Sys, class... Ms>
	StoermerVerlet(Sys&, Ms...)
		-> StoermerVerlet<Sys, monitor::MonitorPack<std::decay_t<Ms>...>>;
}