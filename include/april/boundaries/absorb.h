#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Absorb : Boundary {
		Absorb(): Boundary(-1, false, false) {}

		void apply(env::internal::Particle & particle) const noexcept{
			particle.state = env::ParticleState::DEAD;
		}
	};
}