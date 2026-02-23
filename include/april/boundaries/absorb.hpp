#pragma once

#include "april/boundaries/boundary.hpp"

namespace april {
	struct Absorb : boundary::Boundary {
		static constexpr auto fields = ParticleField::state;

		Absorb(): Boundary(-1, false, false, false) {}

		template<ParticleField M, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> & p, const core::Box &, const Face) const noexcept{
			p.state = ParticleState::DEAD;
		}
	};
}












