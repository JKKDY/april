#pragma once


#include "april/boundaries/boundary.hpp"

namespace april {
	struct Open : boundary::Boundary {
		static constexpr auto fields = ParticleField::none;

		Open(): Boundary(-1, false, false, false) {}

		template<ParticleField M, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> &, const core::Box &, const Face) const noexcept{}
	};
}












