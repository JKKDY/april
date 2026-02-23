#pragma once


#include "april/boundaries/boundary.hpp"

namespace april {
	struct Periodic : boundary::Boundary {
		static constexpr auto fields = ParticleField::position;

		Periodic(): Boundary(-1, true, true, true) {}

		template<ParticleField IncomingMask, particle::IsParticleAttributes U>
		void apply(particle::internal::ScalarParticleRef<IncomingMask, U> & particle, const core::Box & domain_box, const Face face) const noexcept{
			const int sign = boundary::face_sign_pos(face) ? -1 : +1;
			const int ax = boundary::axis_of_face(face);

			particle.position[ax] = particle.position[ax] + sign * domain_box.extent[ax];
		}
	};
}












