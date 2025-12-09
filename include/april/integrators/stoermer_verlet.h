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
			sys.update_all_components();

			// position update
			sys.template for_each_particle<pos_upd_fields>([&](auto p) {
				p.old_position = p.position;
				p.position += dt * p.velocity + (dt*dt) / (2 * p.mass) * p.force;
			}, State::MOVABLE);

			sys.rebuild_structure();
			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			// velocity update
			sys.template for_each_particle<vel_upd_fields>([&](auto p) {
				p.velocity += dt / 2 / p.mass * (p.force + p.old_force);
			}, State::MOVABLE);

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