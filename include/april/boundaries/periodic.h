#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Periodic : Boundary {
		static constexpr env::FieldMask fields = to_field_mask(env::Field::position);

		Periodic(): Boundary(-1, true, true, true) {}

		template<env::IsMutableFetcher F>
	    void apply(F && particle, const env::Box & domain_box, const Face face) const noexcept{
			const int sign = face_sign_pos(face) ? -1 : +1;
			const int ax = axis_of_face(face);

			particle.position()[ax] = particle.position()[ax] + sign * domain_box.extent[ax];
		}
	};
}