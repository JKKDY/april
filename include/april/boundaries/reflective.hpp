#pragma once


#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Reflective : Boundary {
		static constexpr ParticleField fields = ParticleField::position | ParticleField::old_position | ParticleField::velocity;

		Reflective(): Boundary(-1, false, false, true) {}

		template<ParticleField M, core::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<M, U> & particle, const core::Box & domain_box, const Face face) const noexcept{
			const int is_plus = face_sign_pos(face);
			const int ax = axis_of_face(face);

			const vec3 diff = particle.position - particle.old_position;
			vec3 diff_reflected = diff;
			diff_reflected[ax] = -1 * diff_reflected[ax];

			const double y = (is_plus? domain_box.max : domain_box.min)[ax];

			const double t = (y - particle.old_position[ax]) / diff[ax];

			particle.position = particle.old_position + t * diff + (1-t) * diff_reflected;
			particle.velocity[ax] = - particle.velocity[ax];

			AP_ASSERT(particle.position[ax] >= domain_box.min[ax] && particle.position[ax] <= domain_box.max[ax],
				"particle outside of domain on reflected axis! \n\tface:"  + std::to_string(face_to_int(face)) +
				"\n\t" + particle.position.to_string() + "  old pos: " + particle.old_position.to_string() );
		}
	};
}








