#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr FieldMask fields = to_field_mask(Field::none);

		Open(): Boundary(-1, false, false, false) {}

		template<IsUserData UserData>
		void apply(ParticleRef<fields, UserData> &, const Box &, const Face) const noexcept {}	};
}