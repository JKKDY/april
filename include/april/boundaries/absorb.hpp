#pragma once

#include "april/boundaries/boundary.hpp"

namespace april {
	struct AbsorbingBoundary : boundary::Boundary {
		static constexpr auto fields = ParticleField::state;

		AbsorbingBoundary(): Boundary(-1, false, false, false) {}

		template<ParticleField M, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> & p, const core::Box &, const DomainFace) const noexcept{
			p.state = ParticleState::DEAD;
		}
	};
}












