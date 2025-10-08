#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Reflective : Boundary {
		Reflective(): Boundary(-1, false, false) {}

		void apply(env::impl::Particle &) const noexcept{

		}
	};
}