#pragma once

#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Absorb : Boundary {
		static constexpr ParticleField fields = ParticleField::state;

		Absorb(): Boundary(-1, false, false, false) {}

		template<ParticleField M, core::IsUserData U>
		void apply(core::ScalarParticleRef<M, U> & p, const core::Box &, const Face) const noexcept{
			p.state = ParticleState::DEAD;
		}
	};
}




