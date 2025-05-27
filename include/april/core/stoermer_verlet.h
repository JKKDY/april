#pragma once
#include "integrator.h"

namespace april::core {
	class StoermerVerlet : public impl::Integrator<StoermerVerlet> {
		void integration_step() const {
			for (auto &p : env.particles()) {
				p.update_position(dt * p.velocity + pow(dt, 2) / (2 * p.mass) * p.force);
				p.reset_force();
			}

			env.update_forces();

			for (auto &p: env.particles()) {
				p.update_velocity(dt / 2 / p.mass * (p.force + p.old_force));
			}
		}
	};
}