#pragma once


#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr ParticleField fields = ParticleField::none;

		Open(): Boundary(-1, false, false, false) {}

		template<ParticleField M, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> &, const core::Box &, const Face) const noexcept{}
	};
}












