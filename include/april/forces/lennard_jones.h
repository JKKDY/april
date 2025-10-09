#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {

	// Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
	struct LennardJones {
		double epsilon; // Depth of the potential well
		double sigma; // Distance at which potential is zero
		double cutoff_radius; // Maximum interaction distance; negative: no cutoff

		LennardJones(const double epsilon_, const double sigma_, const double cutoff = -1.0)
		: epsilon(epsilon_), sigma(sigma_), sigma2(sigma * sigma) {
			cutoff_radius = (cutoff < 0.0) ? 3.0 * sigma : cutoff;
		}

		vec3 operator()(env::internal::Particle const&, env::internal::Particle const&, vec3 const& r) const noexcept {
			const double r2 = r.norm_squared();
			if (cutoff_radius > 0.0 && r2 > cutoff_radius * cutoff_radius)
				return vec3{0.0, 0.0, 0.0};

			const double inv_r2 = 1.0 / r2;
			const double sigma_r2 = sigma2 * inv_r2;
			const double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
			const double sigma_r12 = sigma_r6 * sigma_r6;
			const double magnitude = 24.0 * epsilon * inv_r2 * (2.0 * sigma_r12 - sigma_r6);

			// Force vector pointing along -r
			return -magnitude * r;
		}

		[[nodiscard]] LennardJones mix(LennardJones const& other) const noexcept {
			// Lorentz-Berthelot mixing rules
			const double mixed_epsilon = std::sqrt(epsilon * other.epsilon);
			const double mixed_sigma = 0.5 * (sigma + other.sigma);
			const double mixed_cutoff = std::sqrt(cutoff_radius * other.cutoff_radius);
			return {mixed_epsilon, mixed_sigma, mixed_cutoff};
		}
	private:
		double sigma2;
	};
}