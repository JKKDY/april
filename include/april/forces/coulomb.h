#pragma once

#include <cmath>

#include "april/common.h"
#include "april/particle/fields.h"
#include "april/forces/force.h"

namespace april::force {
	struct Coulomb : Force{
		double coulomb_constant;

		explicit Coulomb(const double coulomb_const = 1.0, const double cutoff = no_cutoff)
			: Force(cutoff), coulomb_constant(coulomb_const) {}

		template<env::IsConstFetcher F>
		requires requires {
			{ F::user_data_t::charge } -> std::convertible_to<double>;
		}
		vec3 eval(const F & p1, const F & p2, const vec3& r) const noexcept {
			const double inv_r = 1.0 / r.norm();
			const double mag = coulomb_constant * p1.user_data().charge * p2.user_data().charge * inv_r * inv_r;

			return mag * inv_r * r;  // Force vector pointing along +r
		}

		[[nodiscard]] Coulomb mix(Coulomb const& other) const noexcept {
			// Arithmetic average of pre-factor and cutoff
			const double mixed_factor = 0.5 * (coulomb_constant + other.coulomb_constant);
			const double mixed_cutoff = 0.5 * (cutoff + other.cutoff);
			return Coulomb(static_cast<uint8_t>(std::round(mixed_factor)), mixed_cutoff);
		}
	};
}