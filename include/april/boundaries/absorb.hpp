#pragma once

#include "april/particle/fields.hpp"
#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Absorb : Boundary {
		static constexpr env::Field fields = env::Field::state;

		Absorb(): Boundary(-1, false, false, false) {}

		template<env::Field M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> & p, const env::Box &, const Face) const noexcept{
			p.state = env::ParticleState::DEAD;
		}
	};
}

