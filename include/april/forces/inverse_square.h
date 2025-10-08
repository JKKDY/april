#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
 // Inverse-square law force (e.g., gravity, Coulomb). pre_factor scales G or k.
    struct InverseSquare {
        double pre_factor; // e.g. G or Coulomb constant
        double cutoff_radius; // Maximum interaction distance; negative: no cutoff

        explicit InverseSquare(const double pre_factor_ = 1.0, const double cutoff = -1.0)
            : pre_factor(pre_factor_), cutoff_radius(cutoff) {}

        vec3 operator()(env::impl::Particle const& p1, env::impl::Particle const& p2, vec3 const& r) const noexcept {
            const double r2 = r.norm_squared();
            if (cutoff_radius > 0.0 && r2 > cutoff_radius*cutoff_radius) return {};

            const double inv_r = 1.0 / std::sqrt(r2);
            const double inv_r3 = inv_r * inv_r * inv_r;
            const double mag = pre_factor * p1.mass * p2.mass * inv_r3;

            return mag * r;  // Force vector pointing along +r
        }

        [[nodiscard]] InverseSquare mix(InverseSquare const& other) const noexcept {
            // Arithmetic average of prefactor and cutoff
            const double mixed_prefactor = 0.5 * (pre_factor + other.pre_factor);
            const double mixed_cutoff = 0.5 * (cutoff_radius + other.cutoff_radius);
            return InverseSquare(mixed_prefactor, mixed_cutoff);
        }
    };
}