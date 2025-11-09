#pragma once

#include "april/particle/fields.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Absorb : Boundary {
		Absorb(): Boundary(-1, false, false, false) {}

        template<env::IsMutableFetcher F>
		void apply(F && p, const env::Box &, const Face) const noexcept{
			p.state() = env::ParticleState::DEAD;
		}
	};
}