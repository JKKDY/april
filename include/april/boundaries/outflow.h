#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Outflow : Boundary {
		Outflow(): Boundary(-1, false, false) {}

		void apply(env::internal::Particle &) const noexcept{

		}
	};
}