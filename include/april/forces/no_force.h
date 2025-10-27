#pragma once

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce : Force{
		// Negative cutoff_radius means "no cutoff"
		NoForce(): Force(0) {}

		[[nodiscard]] vec3 eval(env::internal::Particle const&, env::internal::Particle const&, vec3 const&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}
	};
}