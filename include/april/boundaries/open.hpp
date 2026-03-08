#pragma once


#include "april/boundaries/boundary.hpp"

namespace april {
	struct OpenBoundary : boundary::Boundary {
		static constexpr auto fields = ParticleField::none;

		OpenBoundary(): Boundary(-1, false, false, false) {}

		void apply(auto, const core::Box &, const DomainFace) const noexcept{}
	};
}














