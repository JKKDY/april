#pragma once

#include "april/particle/fields.hpp"
#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr env::Field fields = env::Field::none;

		Open(): Boundary(-1, false, false, false) {}

		template<env::Field M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> &, const env::Box &, const Face) const noexcept{}
	};
}

