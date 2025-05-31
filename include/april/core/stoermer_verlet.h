#pragma once
#include "april/core/integrator.h"
#include "april/env/particle.h"

namespace april::core {
	template <io::IsOutputWriter OutputW = io::NullOutput>
	class StoermerVerlet : public impl::Integrator<OutputW> {
	public:
		using State = env::ParticleState;
		using Base = impl::Integrator<OutputW>;
		using Base::env;
		using Base::dt;
		using Base::Base;

		void integration_step() const {
			for (auto &p : env.particles(State::MOVABLE)) {
				p.update_position(dt * p.velocity + (dt*dt) / (2 * p.mass) * p.force);
				p.reset_force();
			}

			env.update_forces();

			for (auto &p: env.particles(State::MOVABLE)) {
				p.update_velocity(dt / 2 / p.mass * (p.force + p.old_force));
			}
		}
	};
}