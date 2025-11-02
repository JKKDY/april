#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Open : Boundary {
		using namespace april::env;
		static constexpr env::FieldMask fields = to_field_mask(env::Field::none);

		Open(): Boundary(-1, false, false, false) {}

		template<env::IsUserData UserData>
		void apply(env::ParticleRef<fields, UserData> &, const env::Box &, const Face) const noexcept {}	};
}