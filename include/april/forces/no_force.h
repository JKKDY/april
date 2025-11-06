#pragma once

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		static constexpr env::FieldMask fields = to_field_mask(env::Field::none);
		// Negative cutoff_radius means "no cutoff"
		NoForce(): Force(0) {}


        template<env::IsConstFetcher F>
		vec3 operator()(const F &, const F &, const vec3&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}
	};
}