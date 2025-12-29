#pragma once

#include "april/particle/fields.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Absorb : Boundary {
		static constexpr env::FieldMask fields = +env::Field::state;

		Absorb(): Boundary(-1, false, false, false) {}

		template<env::FieldMask M, env::IsUserData U>
		void apply(env::ParticleRef<M, U> & p, const env::Box &, const Face) const noexcept{
			p.state = env::ParticleState::DEAD;
		}
	};
}