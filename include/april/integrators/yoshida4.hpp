#pragma once
#include "april/integrators/integrator.hpp"
#include "april/core/system.hpp"

#include "april/monitors/monitor.hpp"

namespace april {

	template <core::IsSystem Sys, monitor::internal::IsMonitorPack Monitors>
	class Yoshida4 : public integrator::Integrator<Sys, Monitors> {
	public:
		using State = ParticleState;
		using Base = integrator::Integrator<Sys, Monitors>;
		using Base::dt;
		using Base::sys;
		using Base::Base;

		static constexpr ParticleField pos_upd_fields =
			ParticleField::state | ParticleField::velocity | ParticleField::position | ParticleField::mass |
			ParticleField::old_position | ParticleField::force;

		static constexpr ParticleField vel_upd_fields =
			ParticleField::state | ParticleField::velocity | ParticleField::force | ParticleField::mass;

		void velocity_verlet_step(double delta_t) const {
			sys.update_all_components();

			sys.template for_each_particle<Sys::parallel_policy>(
				april::universal_kernel<pos_upd_fields, pos_upd_fields>([&](auto p) {
					p.old_position = p.position;
					p.velocity += (delta_t / 2.0) * (p.force / p.mass);
					p.position += delta_t * p.velocity;
				}
			), State::MOVABLE);

			sys.rebuild_structure();
			sys.apply_boundary_conditions();
			sys.update_forces();
			sys.apply_force_fields();

			sys.template for_each_particle<Sys::parallel_policy>(
				april::universal_kernel<vel_upd_fields, vel_upd_fields>([&](auto p) {
					p.velocity += (delta_t / 2.0) * (p.force / p.mass);
				}
			), State::MOVABLE);

			sys.apply_controllers();
		}

		void integration_step() const {
			constexpr double w1 = 1.3512071919596578;
			constexpr double w2 = -1.7024143839193153;

			velocity_verlet_step(w1*dt);
			velocity_verlet_step(w2*dt);
			velocity_verlet_step(w1*dt);
		}
	};

	// Deduction guide so user can write Yoshida4(sys, MonitorPack<M1, M2, M3>)
	template<core::IsSystem Sys, class... Ms>
	Yoshida4(Sys&, monitor::internal::MonitorPack<Ms...>)
		-> Yoshida4<Sys, monitor::internal::MonitorPack<Ms...>>;

	// Deduction guide so user can write StoermerVerlet(sys, m1, m2, m3)
	template<core::IsSystem Sys, class... Ms>
	requires (monitor::IsMonitor<std::decay_t<Ms>> && ...)
	Yoshida4(Sys&, Ms...)
		-> Yoshida4<Sys, monitor::internal::MonitorPack<std::decay_t<Ms>...>>;
}














