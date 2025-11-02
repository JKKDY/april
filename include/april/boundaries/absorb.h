#pragma once

#include "april/env/particle.h"
#include "april/boundaries/boundary.h"

namespace april::boundary {
	struct Absorb : Boundary {
		static constexpr env::FieldMask fields = env::to_field_mask(env::Field::state);
		Absorb(): Boundary(-1, false, false, false) {}

		template<env::IsUserData UserData>
		void apply(env::ParticleRef<fields, UserData> p, const env::Box &, const Face) const noexcept{
			p.state = env::ParticleState::DEAD;
		}
	};
}