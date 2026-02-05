#pragma once

#include "april/base/types.hpp"
#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"


namespace april::force {
	// Harmonic spring force (Hooke's law). k: spring constant; r0: equilibrium distance.
	struct Harmonic : Force{
		static constexpr env::FieldMask fields = +env::Field::none;

		double k; // Spring constant
		double r0; // Equilibrium distance

		Harmonic(const double strength, const double equilibrium, const double cutoff = no_cutoff)
		: Force(cutoff), k(strength), r0(equilibrium) {}


		auto eval(auto, auto, const auto& r) const noexcept {
			const auto dist = r.norm();
			const auto magnitude = k * (dist - r0) / dist; // F = k * (dist - r0) * (r / dist)
			return magnitude * r;
		}

		[[nodiscard]] Harmonic mix(Harmonic const& other) const noexcept {
			// Arithmetic average of k and r0; carry max cutoff
			const double mixed_k = 2* k * other.k / (k + other.k); // models 2 springs in series
			const double mixed_r0 = 0.5 * (r0 + other.r0);

			const Harmonic h(mixed_k, mixed_r0, std::max(cutoff(), other.cutoff()));
			return h;
		}

		bool operator==(const Harmonic&) const = default;
	};
}