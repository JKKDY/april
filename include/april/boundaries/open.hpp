#pragma once

#include "april/particle/fields.hpp"
#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr ParticleField fields = ParticleField::none;

		Open(): Boundary(-1, false, false, false) {}

		template<ParticleField M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> &, const env::Box &, const Face) const noexcept{}
	};
}



