#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// Harmonic spring force (Hooke's law). k: spring constant; r0: equilibrium distance.
	struct Harmonic : Force{
		double k; // Spring constant
		double r0; // Equilibrium distance

		Harmonic(const double strength, const double equilibrium, const double cutoff = -1)
		: Force(cutoff), k(strength), r0(equilibrium) {}

		vec3 operator()(env::internal::Particle const&, env::internal::Particle const&, vec3 const& r) const noexcept {
			const double dist = r.norm();
			if (dist > cutoff) return {};
			const double magnitude = k * (dist - r0) / dist; // F = -k * (dist - r0) * (r / dist)
			return -magnitude * r;
		}

		[[nodiscard]] Harmonic mix(Harmonic const& other) const noexcept {
			// Arithmetic average of k and r0; carry max cutoff
			const double mixed_k = 0.5 * (k + other.k);
			const double mixed_r0 = 0.5 * (r0 + other.r0);
			Harmonic h(mixed_k, mixed_r0);
			h.cutoff = std::max(cutoff, other.cutoff);
			return h;
		}
	};
}