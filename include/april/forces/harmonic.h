#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
	// Harmonic spring force (Hooke's law). k: spring constant; r0: equilibrium distance.
	struct Harmonic {
		double k; // Spring constant
		double r0; // Equilibrium distance
		double cutoff_radius; // Negative: no cutoff

		Harmonic(const double k_, const double r0_)
		: k(k_), r0(r0_), cutoff_radius(-1.0) {}

		vec3 operator()(env::internal::Particle const&, env::internal::Particle const&, vec3 const& r) const noexcept {
			const double dist = r.norm();
			const double magnitude = k * (dist - r0) / dist; // F = -k * (dist - r0) * (r / dist)
			return -magnitude * r;
		}

		[[nodiscard]] Harmonic mix(Harmonic const& other) const noexcept {
			// Arithmetic average of k and r0; carry max cutoff
			const double mixed_k = 0.5 * (k + other.k);
			const double mixed_r0 = 0.5 * (r0 + other.r0);
			Harmonic h(mixed_k, mixed_r0);
			h.cutoff_radius = std::max(cutoff_radius, other.cutoff_radius);
			return h;
		}
	};
}