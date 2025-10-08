#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Periodic : Boundary {
		Periodic(): Boundary(-1, true, true) {}

		void apply(env::impl::Particle &) const noexcept{

		}
	};
}