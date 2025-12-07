#pragma once

#include "april/particle/fields.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr env::FieldMask fields = +env::Field::none;

		Open(): Boundary(-1, false, false, false) {}

		template<env::FieldMask M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> &, const env::Box &, const Face) const noexcept{}
	};
}