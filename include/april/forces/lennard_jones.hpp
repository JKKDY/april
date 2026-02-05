#pragma once

#include <cmath>
#include "april/forces/force.hpp"
#include "april/base/types.hpp"
#include "april/particle/fields.hpp"


namespace april::force {

	// Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
	struct LennardJones : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		LennardJones(const double epsilon, const double sigma, const double cutoff = -1.0)
		: Force(cutoff < 0.0 ? 3.0 * sigma : cutoff), epsilon_(epsilon), sigma_(sigma) {
			calculate_constants();
		}

		auto epsilon(const double e)  {
			epsilon_ = e;
			calculate_constants();
			return *this;
		}

		auto sigma(const double s)  {
			sigma_ = s;
			calculate_constants();
			return *this;
		}

		auto eval(auto, auto, const auto& r) const noexcept {
			const auto inv_r2 = static_cast<vec3::type>(1.0) / (r.x*r.x + r.y*r.y + r.z * r.z);
			const auto inv_r6 = inv_r2 * inv_r2 * inv_r2;
			const auto magnitude = (c12_force * inv_r6 - c6_force) * inv_r6 * inv_r2;
			//
			// // Force vector pointing along -r
			return -magnitude * r;
		}

		[[nodiscard]] LennardJones mix(LennardJones const& other) const noexcept {
			// Lorentz-Berthelot mixing rules
			const double mixed_epsilon = std::sqrt(epsilon_ * other.epsilon_);
			const double mixed_sigma = 0.5 * (sigma_ + other.sigma_);
			const double mixed_cutoff = std::sqrt(cutoff() * other.cutoff());
			return {mixed_epsilon, mixed_sigma, mixed_cutoff};
		}

		bool operator==(const LennardJones&) const = default;

	private:
		void calculate_constants() {
			const vec3::type sigma2 = sigma_ * sigma_;
			const vec3::type sigma6 = sigma2 * sigma2 * sigma2;
			const vec3::type sigma12 = sigma6 * sigma6;

			c6_force = 24.0 * epsilon_ * sigma6;
			c12_force = 48.0 * epsilon_ * sigma12;
		}

		// Precomputed force constants
		vec3::type c12_force;
		vec3::type c6_force;

		double epsilon_; // Depth of the potential well
		double sigma_; // Distance at which potential is zero
	};
}