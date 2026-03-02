#pragma once

#include <cmath>

#include "april/base/types.hpp"
#include "april/forces/force.hpp"

namespace april {
	struct Coulomb : force::Force{
		static constexpr auto fields = ParticleField::attributes;

		double coulomb_constant;

		explicit Coulomb(const double coulomb_const = 1.0, const double cutoff = force::no_cutoff)
			: Force(cutoff), coulomb_constant(coulomb_const) {}

		template<ParticleField M, ParticleField N, particle::IsParticleAttributes U>
		requires requires { // TODO replace with concept that checks if particle has U::charge
			{ U::charge } -> std::convertible_to<double>;
		}
		auto eval(const particle::internal::ScalarParticleView<M, N, U> & p1, const particle::internal::ScalarParticleView<M, N, U> & p2, const vec3& r) const noexcept {
			const double inv_r = r.inv_norm();
			const double mag = coulomb_constant * p1.attributes.charge * p2.attributes.charge * inv_r * inv_r;

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












