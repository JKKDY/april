#pragma once

#include <cmath>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
	#include <immintrin.h>
#endif
#include "april/common.hpp"
#include "april/particle/fields.hpp"


namespace april::force {

	// Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
	struct LennardJones : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		static double fast_inv_r2(const double r_squared) {
#if defined(__x86_64__) || defined(_M_X64)
			// Intel optimized path (Approximation with relative error <= 2^-14)
			const __m128d val = _mm_set_sd(r_squared);
			const __m128d result = _mm_rcp14_sd(val, val);
			return _mm_cvtsd_f64(result);
#else
			// ARM/Apple Silicon & Generic Fallback
			// The FPU division on Apple M-series chips is extremely fast
			// often making the approximation unnecessary.
			return 1.0 / r_squared;
#endif
		}
		LennardJones(const double epsilon_, const double sigma_, const double cutoff = -1.0)
		: Force(cutoff < 0.0 ? 3.0 * sigma_ : cutoff), epsilon(epsilon_), sigma(sigma_) {
			const vec3::type sigma2 = sigma * sigma;
			const vec3::type sigma6 = sigma2 * sigma2 * sigma2;
			const vec3::type sigma12 = sigma6 * sigma6;

			c6_force = 24.0 * epsilon * sigma6;
			c12_force = 48.0 * epsilon * sigma12;
		}

		template<env::FieldMask M, env::IsUserData U>
		vec3 eval(const env::ParticleView<M, U> &, const env::ParticleView<M, U> &, const vec3& r) const noexcept {

			// const vec3::type inv_r2 = static_cast<vec3::type>(1.0) / (r.x*r.x + r.y*r.y + r.z * r.z);
			const vec3::type inv_r2 = fast_inv_r2(r.norm_squared());
			const vec3::type inv_r6 = inv_r2 * inv_r2 * inv_r2;
			const vec3::type magnitude = (c12_force * inv_r6 - c6_force) * inv_r6 * inv_r2;

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
		vec3::type c12_force;
		vec3::type c6_force;

		double epsilon; // Depth of the potential well
		double sigma; // Distance at which potential is zero
	};
}