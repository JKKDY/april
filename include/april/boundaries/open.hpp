#pragma once


#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr ParticleField fields = ParticleField::none;

		Open(): Boundary(-1, false, false, false) {}

		template<ParticleField M, core::IsUserData U>
		void apply(core::ScalarParticleRef<M, U> &, const core::Box &, const Face) const noexcept{}
	};
}




