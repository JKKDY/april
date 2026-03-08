#pragma once

#include "april/boundaries/boundary.hpp"

namespace april {
	struct AbsorbingBoundary : boundary::Boundary {
		static constexpr auto fields = ParticleField::state;

		AbsorbingBoundary(): Boundary(-1, false, false, false) {}

		void apply(auto && p, const core::Box &, const DomainFace) const noexcept{
			p.state = ParticleState::DEAD;
		}
	};
}














