#pragma once


#include "april/boundaries/boundary.hpp"

namespace april {
	struct PeriodicBoundary : boundary::Boundary {
		static constexpr auto fields = ParticleField::position;

		PeriodicBoundary(): Boundary(-1, true, true, true) {}

		void apply(auto && particle, const core::Box & domain_box, const DomainFace face) const noexcept{
			const int sign = boundary::face_sign_pos(face) ? -1 : +1;
			const int ax = boundary::axis_of_face(face);

			particle.position[ax] = particle.position[ax] + sign * domain_box.extent[ax];
		}
	};
}












