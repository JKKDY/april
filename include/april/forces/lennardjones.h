#pragma once

#include <cmath>

#include "april/common.h"
#include "april/particle/fields.h"


namespace april::force {

	// Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
	struct LennardJones : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		double epsilon; // Depth of the potential well
		double sigma; // Distance at which potential is zero

		LennardJones(const double epsilon_, const double sigma_, const double cutoff = -1.0)
		: Force(cutoff < 0.0 ? 3.0 * sigma_ : cutoff),
		epsilon(epsilon_), sigma(sigma_), sigma2(sigma * sigma) {}


		template<env::IsUserData U>
		vec3 eval(const env::ParticleView<fields, U> &, const env::ParticleView<fields, U> &, const vec3& r) const noexcept {
			const double inv_r2 = 1.0 / r.norm_squared();
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
			const double mixed_cutoff = std::sqrt(cutoff() * other.cutoff());
			return {mixed_epsilon, mixed_sigma, mixed_cutoff};
		}
	private:
		double sigma2;
	};
}