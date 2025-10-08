#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	template <class Force>
	struct Repulsive : Boundary {
		explicit Repulsive(Force & force): Boundary(force.cutoff(), false, false) {}

		void apply(env::impl::Particle &) const noexcept{

		}
	};
}