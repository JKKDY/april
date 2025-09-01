#pragma once
#include "april/core/integrator.h"
#include "april/core/system.h"
#include "april/env/particle.h"
#include "april/io/monitor.h"
#include "april/defaults.h"

namespace april::core {

	template<IsSystem Sys, class Pack> class Yoshida4;

	template <IsSystem Sys, class ... TMonitors>
	class Yoshida4<Sys, io::MonitorPack<TMonitors...>>
		: public impl::Integrator<Sys, io::MonitorPack<TMonitors...>> {
	public:
		using State = env::ParticleState;
		using Base = impl::Integrator<Sys, io::MonitorPack<TMonitors...>>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		void stoermer_verlet_step(double delta_t) const {
			for (auto i = sys.index_start(); i <= sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_position(delta_t * p.velocity + (delta_t*delta_t) / (2 * p.mass) * p.force);
			}

			sys.update_forces();

			for (auto i = sys.index_start(); i <= sys.index_end(); ++i) {
				auto & p = sys.get_particle_by_index(i);
				if (static_cast<int>(p.state & State::MOVABLE))
					p.update_velocity(delta_t / 2 / p.mass * (p.force + p.old_force));
			}
		}

		void integration_step() const {
			constexpr double w1 = 1.3512071919596578;
			constexpr double w2 = -1.7024143839193153;

			stoermer_verlet_step(w1*dt);
			stoermer_verlet_step(w2*dt);
			stoermer_verlet_step(w1*dt);
		}
	};

	// Deduction guide so user can write StoermerVerlet(sys, MonitorPack<M1, M2, M3>)
	template<class Sys, class... Ms>
	Yoshida4(Sys&, io::MonitorPack<Ms...>)
		-> Yoshida4<Sys, io::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys)
	template<class Sys>
	Yoshida4(Sys&)
		-> Yoshida4<Sys, DefaultMonitors>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<class Sys, class... Ms>
	Yoshida4(Sys&, Ms...)
		-> Yoshida4<Sys, io::MonitorPack<std::decay_t<Ms>...>>;
}