#pragma once
#include "april/core/integrator.h"
#include "april/env/particle.h"

namespace april::core {
	class StoermerVerlet : public impl::Integrator<> {
		using State = env::ParticleState;
		void integration_step() const {
			for (auto &p : env.particles(State::MOVABLE)) {
				p.update_position(dt * p.velocity + pow(dt, 2) / (2 * p.mass) * p.force);
				p.reset_force();
			}

			env.update_forces();

			for (auto &p: env.particles(State::MOVABLE)) {
				p.update_velocity(dt / 2 / p.mass * (p.force + p.old_force));
			}
		}
	};
}