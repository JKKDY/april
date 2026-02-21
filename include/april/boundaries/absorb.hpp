#pragma once

#include "april/particle/fields.hpp"
#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Absorb : Boundary {
		static constexpr ParticleField fields = ParticleField::state;

		Absorb(): Boundary(-1, false, false, false) {}

		template<ParticleField M, env::IsUserData U>
		void apply(env::ScalarParticleRef<M, U> & p, const env::Box &, const Face) const noexcept{
			p.state = ParticleState::DEAD;
		}
	};
}



