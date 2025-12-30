#pragma once
#include "april/integrators/integrator.hpp"
#include "april/system/system.hpp"
#include "april/particle/fields.hpp"
#include "april/monitors/monitor.hpp"

namespace april::integrator {

	template<core::IsSystem Sys, class Pack> class VelocityVerlet;

	template <core::IsSystem Sys, class ... TMonitors>
	class VelocityVerlet<Sys, monitor::MonitorPack<TMonitors...>>
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
			env::Field::state | env::Field::velocity | env::Field::force | env::Field::mass;

		void integration_step() const {
			sys.update_all_components();

			sys.template for_each_particle<pos_upd_fields>([&](auto p) {
				p.old_position = p.position;
				p.velocity += (dt / 2.0) * (p.force / p.mass);
				p.position += dt * p.velocity;
			}, State::MOVABLE);

			sys.rebuild_structure();
			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			sys.template for_each_particle<vel_upd_fields>([&](auto p) {
			   p.velocity += (dt / 2.0) * (p.force / p.mass);
			}, State::MOVABLE);

			sys.apply_controllers();
		}
	};

	// Deduction guide so user can write StoermerVerlet(sys, MonitorPack<M1, M2, M3>)
	template<class Sys, class... Ms>
	VelocityVerlet(Sys&, monitor::MonitorPack<Ms...>)
		-> VelocityVerlet<Sys, monitor::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<class Sys, class... Ms>
	VelocityVerlet(Sys&, Ms...)
		-> VelocityVerlet<Sys, monitor::MonitorPack<std::decay_t<Ms>...>>;
}