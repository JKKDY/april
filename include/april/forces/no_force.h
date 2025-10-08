#pragma once

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// No-op force: always returns zero vector and mixes to itself.
	struct NoForce {
		// Negative cutoff_radius means "no cutoff"
		double cutoff_radius = 0.0;

		vec3 operator()(env::impl::Particle const&, env::impl::Particle const&, vec3 const&) const noexcept {
			return vec3{0.0, 0.0, 0.0};
		}

		[[nodiscard]] NoForce mix(NoForce const&) const noexcept {
			return {};
		}
	};
}