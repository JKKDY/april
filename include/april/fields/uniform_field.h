#pragma once

#include "april/fields/field.h"
#include "april/env/particle.h"

namespace april::field {
	struct UniformField final : Field {

		UniformField(const vec3 & force_dir): force(force_dir) {}

		void apply(env::RestrictedParticleRef particle) const {
			particle.force += force;
		}

	private:
		const vec3 force;
	};
}