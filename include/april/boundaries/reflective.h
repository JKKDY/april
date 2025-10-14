#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Reflective : Boundary {
		Reflective(): Boundary(-1, false, false, true) {}

		void apply(env::internal::Particle & particle, const env::Box & domain_box, const Face face) const noexcept{
			const int is_plus = face_sign_pos(face);
			const int ax = axis_of_face(face);

			// AP_ASSERT(domain_box.contains(particle.old_position),
			//           "Domain does not contain particle!: old_pos" + particle.old_position.to_string());


			// if (particle.id == 1497) {
			// 	std::cout << "reflecting particle 1497\n";
			// 	std::cout << "pos:     " + particle.position.to_string() + "\n";
			// 	std::cout << "oldpos:  " + particle.old_position.to_string() + "\n";
			// 	std::cout << "vel:     " + particle.velocity.to_string() + "\n";
			//
			// 	const vec3 diff = particle.position - particle.old_position;
			// 	vec3 diff_reflected = diff;
			// 	diff_reflected[ax] = -1 * diff_reflected[ax];
			//
			// 	std::cout << "diff:    " + diff.to_string() + "\n";
			// 	std::cout << "rdiff:   " + diff_reflected.to_string() + "\n";
			//
			// 	const double y = (is_plus? domain_box.max : domain_box.min)[ax];
			// 	const double t = (y - particle.old_position[ax]) / diff[ax];
			//
			// 	std::cout << "y:       " + std::to_string(y) + "\n";
			// 	std::cout << "t:       " + std::to_string(t) + "\n";
			//
			// 	const vec3 new_pos = particle.old_position + t * diff + (1-t) * diff_reflected;
			// 	vec3 new_vel = particle.velocity;
			// 	new_vel[ax] = - particle.velocity[ax];
			//
			// 	std::cout << "new_pos: " + new_pos.to_string() + "\n";
			// 	std::cout << "new_vel: " + new_vel.to_string() + "\n";
			// }


			const vec3 diff = particle.position - particle.old_position;
			vec3 diff_reflected = diff;
			diff_reflected[ax] = -1 * diff_reflected[ax];

			const double y = (is_plus? domain_box.max : domain_box.min)[ax];

			const double t = (y - particle.old_position[ax]) / diff[ax];
			// AP_ASSERT(t >= 0 && t < 1, "t should be between 0 and 1, got" + std::to_string(t) +
			// 	"on face: " + std::to_string(face_to_int(face)) +
			// 	"\n\t pos: " + particle.position.to_string() + "  old pos: " + particle.old_position.to_string() );

			particle.position = particle.old_position + t * diff + (1-t) * diff_reflected;
			particle.velocity[ax] = - particle.velocity[ax];

			AP_ASSERT(particle.position[ax] >= domain_box.min[ax] && particle.position[ax] <= domain_box.max[ax],
				"particle outside of domain on reflected axis! \n\tface:"  + std::to_string(face_to_int(face)) +
				"\n\t" + particle.position.to_string() + "  old pos: " + particle.old_position.to_string() );
		}
	};
}