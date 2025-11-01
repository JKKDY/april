#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Periodic : Boundary {
		static constexpr FieldMask fields = to_field_mask(Field::position);

		Periodic(): Boundary(-1, true, true, true) {}

		template<IsUserData UserData>
		void apply(ParticleRef<fields, UserData> particle, const Box & domain_box, const Face face) const noexcept{
			const int sign = face_sign_pos(face) ? -1 : +1;
			const int ax = axis_of_face(face);

			particle.position[ax] = particle.position[ax] + sign * domain_box.extent[ax];
		}
	};
}