#pragma once

#include "april/particle/fields.hpp"
#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Periodic : Boundary {
		static constexpr env::FieldMask fields = +env::Field::position;

		Periodic(): Boundary(-1, true, true, true) {}

		template<env::FieldMask IncomingMask, env::IsUserData U>
		void apply(env::ParticleRef<IncomingMask, U> & particle, const env::Box & domain_box, const Face face) const noexcept{
			const int sign = face_sign_pos(face) ? -1 : +1;
			const int ax = axis_of_face(face);

			particle.position[ax] = particle.position[ax] + sign * domain_box.extent[ax];
		}
	};
}