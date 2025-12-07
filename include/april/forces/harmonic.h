#pragma once

#include <cmath>

#include "april/common.h"
#include "april/particle/fields.h"
#include "april/forces/force.h"



namespace april::force {
	// Harmonic spring force (Hooke's law). k: spring constant; r0: equilibrium distance.
	struct Harmonic : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		double k; // Spring constant
		double r0; // Equilibrium distance

		Harmonic(const double strength, const double equilibrium, const double cutoff = no_cutoff)
		: Force(cutoff), k(strength), r0(equilibrium) {}


		template<env::FieldMask M, env::IsUserData U>
		vec3 eval(const env::ParticleView<M, U> &, const env::ParticleView<M, U> &, const vec3& r) const noexcept {
			const double dist = r.norm();
			const double magnitude = k * (dist - r0) / dist; // F = -k * (dist - r0) * (r / dist)
			return -magnitude * r;
		}

		[[nodiscard]] Harmonic mix(Harmonic const& other) const noexcept {
			// Arithmetic average of k and r0; carry max cutoff
			const double mixed_k = 0.5 * (k + other.k);
			const double mixed_r0 = 0.5 * (r0 + other.r0);
			const Harmonic h(mixed_k, mixed_r0, std::max(cutoff(), other.cutoff()));
			return h;
		}
	};
}