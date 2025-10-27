#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Open : Boundary {
		Open(): Boundary(-1, false, false, false) {}

		void apply(env::internal::Particle &, const env::Box &, const Face) const noexcept {}
	};
}