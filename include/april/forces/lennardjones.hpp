#pragma once

#include <cmath>

#include "april/common.hpp"
#include "april/particle/fields.hpp"


namespace april::force {

	// Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
	struct LennardJones : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		LennardJones(const double epsilon_, const double sigma_, const double cutoff = -1.0)
		: Force(cutoff < 0.0 ? 3.0 * sigma_ : cutoff), epsilon(epsilon_), sigma(sigma_) {
			const double sigma2 = sigma * sigma;
			const double sigma6 = sigma2 * sigma2 * sigma2;
			const double sigma12 = sigma6 * sigma6;

			c6_force = 24.0 * epsilon * sigma6;
			c12_force = 48.0 * epsilon * sigma12;
		}

		template<env::FieldMask M, env::IsUserData U>
		vec3 eval(const env::ParticleView<M, U> &, const env::ParticleView<M, U> &, const vec3& r) const noexcept {
			const double inv_r = r.inv_norm();
			const double inv_r2 = inv_r * inv_r;
			const double inv_r6 = inv_r2 * inv_r2 * inv_r2;

			const double magnitude = (c12_force * inv_r6 - c6_force) * inv_r6 * inv_r2;

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

		bool operator==(const LennardJones&) const = default;
	private:
		// Precomputed force constants
		double c12_force;
		double c6_force;

		double epsilon; // Depth of the potential well
		double sigma; // Distance at which potential is zero
	};
}