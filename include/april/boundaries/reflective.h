#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Reflective : Boundary {
		Reflective(): Boundary(-1, false, false, true) {}

		void apply(env::internal::Particle & particle, const env::Box & domain_box, const Face face) const noexcept{
			const int is_plus = face_sign_pos(face);
			const int ax = axis_of_face(face);

			const vec3 diff = particle.position - particle.old_position;
			vec3 diff_reflected = diff;
			diff_reflected[ax] = -1 * diff_reflected[ax];

			const double y = (is_plus? domain_box.max : domain_box.min)[ax];

			const double t = (y - particle.position[ax]) / diff[ax];
			AP_ASSERT(t >= 0 && t < 1, "t should be between 0 and 1");

			particle.position = particle.old_position + t * diff + (1-t) * diff_reflected;
			particle.velocity[ax] = - particle.velocity[ax];
		}
	};
}