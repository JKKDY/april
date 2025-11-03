#pragma once

#include "april/fields/field.h"
#include "april/env/particle.h"

namespace april::field {
	struct UniformField final : Field {
		static constexpr env::FieldMask fields = to_field_mask(env::Field::force);

		explicit UniformField(const vec3 & force_dir): force(force_dir) {}

		template<env::IsUserData U>
		void apply(env::RestrictedParticleRef<fields, U> particle) const {
			particle.force += force;
		}

	private:
		const vec3 force;
	};
}