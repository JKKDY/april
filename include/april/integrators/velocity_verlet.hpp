#pragma once
#include "april/integrators/integrator.hpp"
#include "april/core/system.hpp"

#include "april/monitors/monitor.hpp"

namespace april::integrator {

	template<core::IsSystem Sys, monitor::internal::IsMonitorPack Pack> class VelocityVerlet;

	template <core::IsSystem Sys,  monitor::internal::IsMonitorPack Monitors>
	class VelocityVerlet
		: public Integrator<Sys, Monitors> {
	public:
		using State = ParticleState;
		using Base = Integrator<Sys, Monitors>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		static constexpr ParticleField pos_upd_fields =
			ParticleField::state | ParticleField::velocity | ParticleField::position | ParticleField::mass |
			ParticleField::old_position | ParticleField::force;

		static constexpr ParticleField vel_upd_fields =
			ParticleField::state | ParticleField::velocity | ParticleField::force | ParticleField::mass;

		void integration_step() const {
			sys.update_all_components();

			sys.template for_each_particle<pos_upd_fields>(april::universal_kernel(
				[&](auto p) {
					p.old_position = p.position;
					p.velocity += (dt / 2.0) * (p.force / p.mass);
					p.position += dt * p.velocity;
				}
			), State::MOVABLE);

			sys.rebuild_structure();
			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			sys.template for_each_particle<vel_upd_fields>(april::universal_kernel(
				[&](auto p) {
					p.velocity += (dt / 2.0) * (p.force / p.mass);
				}
			), State::MOVABLE);

			sys.apply_controllers();
		}
	};

	// Deduction guide so user can write StoermerVerlet(sys, MonitorPack<M1, M2, M3>)
	template<core::IsSystem Sys, monitor::IsMonitor... Ms>
	VelocityVerlet(Sys&, monitor::internal::MonitorPack<Ms...>)
		-> VelocityVerlet<Sys, monitor::internal::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<core::IsSystem Sys,typename ... Ms>
	requires (monitor::IsMonitor<std::decay_t<Ms>> && ...)
	VelocityVerlet(Sys&, Ms...)
		-> VelocityVerlet<Sys, monitor::internal::MonitorPack<std::decay_t<Ms>...>>;
}

