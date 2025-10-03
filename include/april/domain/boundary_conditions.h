#pragma once

#include "april/env/boundary.h"
#inlcude "april/env/particle.h"

namespace april::env {
	class Outflow : public BoundaryCondition {
		void apply(env::impl::Particle &) const noexcept{
			// No-Op
		}
	};
}


