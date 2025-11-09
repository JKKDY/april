#pragma once

#include "april/particle/particle_fields.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Open : Boundary {
		static constexpr env::FieldMask fields = to_field_mask(env::Field::none);

		Open(): Boundary(-1, false, false, false) {}

		template<env::IsMutableFetcher F>
	    void apply(F &&, const env::Box &, const Face) const noexcept {}
	};
}