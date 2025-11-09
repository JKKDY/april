#pragma once
#include "april/integrators/integrator.h"
#include "april/system/system.h"
#include "april/particle/fields.h"
#include "april/monitors/monitor.h"

namespace april::integrator {

	template<core::IsSystem Sys, class Pack> class StoermerVerlet;

	template <core::IsSystem Sys, class ... TMonitors>
	class StoermerVerlet<Sys, monitor::MonitorPack<TMonitors...>>
		: public Integrator<Sys, monitor::MonitorPack<TMonitors...>> {
	public:
		using State = env::ParticleState;
		using Base = Integrator<Sys, monitor::MonitorPack<TMonitors...>>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		static constexpr env::FieldMask pos_upd_fields =
			env::Field::state | env::Field::velocity | env::Field::position | env::Field::mass | env::Field::old_position | env::Field::force;

		static constexpr env::FieldMask vel_upd_fields =
			env::Field::state | env::Field::velocity | env::Field::force | env::Field::mass | env::Field::old_force;

		void integration_step() const {
			// position update
			for (auto i = sys.index_start(); i < sys.index_end(); ++i) {
				auto p = sys.template get_particle_by_index<pos_upd_fields>(i);
				if (static_cast<int>(p.state & State::MOVABLE)) {
					p.old_position = p.position;
					p.position += dt * p.velocity + (dt*dt) / (2 * p.mass) * p.force;
				}
			}

			sys.register_all_particle_movements();
			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			// velocity update
			for (auto i = sys.index_start(); i < sys.index_end(); ++i) {
				auto p = sys.template get_particle_by_index<vel_upd_fields>(i);
				if (static_cast<int>(p.state & State::MOVABLE)) {
					p.velocity += dt / 2 / p.mass * (p.force + p.old_force);
				}
			}

			sys.apply_controllers();
		}
	};

	// Deduction guide so user can write StoermerVerlet(sys, MonitorPack<M1, M2, M3>)
	template<class Sys, class... Ms>
	StoermerVerlet(Sys&, monitor::MonitorPack<Ms...>)
		-> StoermerVerlet<Sys, monitor::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<class Sys, class... Ms>
	StoermerVerlet(Sys&, Ms...)
		-> StoermerVerlet<Sys, monitor::MonitorPack<std::decay_t<Ms>...>>;
}