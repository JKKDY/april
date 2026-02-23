#pragma once


#include "april/boundaries/boundary.hpp"

namespace april {
	struct OpenBoundary : boundary::Boundary {
		static constexpr auto fields = ParticleField::none;

		OpenBoundary(): Boundary(-1, false, false, false) {}

		template<ParticleField M, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> &, const core::Box &, const DomainFace) const noexcept{}
	};
}












