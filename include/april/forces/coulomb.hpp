#pragma once

#include <cmath>

#include "april/base/types.hpp"
#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"

namespace april::force {
	struct Coulomb : Force{
		static constexpr env::FieldMask fields = +env::Field::user_data;

		double coulomb_constant;

		explicit Coulomb(const double coulomb_const = 1.0, const double cutoff = no_cutoff)
			: Force(cutoff), coulomb_constant(coulomb_const) {}

		template<env::FieldMask M, env::IsUserData U>
		requires requires {
			{ U::charge } -> std::convertible_to<double>;
		}
		auto eval(const env::ParticleView<M, U> & p1, const env::ParticleView<M, U> & p2, const vec3& r) const noexcept {
			const double inv_r = r.inv_norm();
			const double mag = coulomb_constant * p1.user_data.charge * p2.user_data.charge * inv_r * inv_r;

			return mag * inv_r * r;  // Force vector pointing along +r
		}

		[[nodiscard]] Coulomb mix(Coulomb const& other) const {
			// Arithmetic average of pre-factor and cutoff
			if (std::abs(coulomb_constant - other.coulomb_constant) > 1e-9) {
				throw std::invalid_argument("Cannot mix different Coulomb Constants!");
			}
			return Coulomb(coulomb_constant, std::max(cutoff(), other.cutoff()));
		}

		bool operator==(const Coulomb&) const = default;
	};

}