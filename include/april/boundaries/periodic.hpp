#pragma once


#include "april/boundaries/boundary.hpp"

namespace april::boundary {
	struct Periodic : Boundary {
		static constexpr ParticleField fields = ParticleField::position;

		Periodic(): Boundary(-1, true, true, true) {}

		template<ParticleField IncomingMask, core::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<IncomingMask, U> & particle, const core::Box & domain_box, const Face face) const noexcept{
			const int sign = face_sign_pos(face) ? -1 : +1;
			const int ax = axis_of_face(face);

			particle.position[ax] = particle.position[ax] + sign * domain_box.extent[ax];
		}
	};
}








