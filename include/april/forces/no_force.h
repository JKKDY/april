#pragma once

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		static constexpr env::FieldMask fields = env::to_field_mask(env::Field::none);
		// Negative cutoff_radius means "no cutoff"
		NoForce(): Force(0) {}


		template<env::IsUserData U1, env::IsUserData U2>
		vec3 operator()(env::ParticleView<fields, U1>, env::ParticleView<fields, U2>, const vec3&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}
	};
}