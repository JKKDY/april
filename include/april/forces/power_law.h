#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
    // Power force law force (e.g., gravity, Coulomb)

    struct PowerLaw {
        double pre_factor; // e.g. G or Coulomb constant
        double cutoff_radius; // Maximum interaction distance; negative: no cutoff
        uint8_t exp;

        explicit PowerLaw(const uint8_t exp, const double pre_factor_ = 1.0, const double cutoff = -1.0)
            : pre_factor(pre_factor_), cutoff_radius(cutoff), exp(exp) {}

        vec3 operator()(env::internal::Particle const& p1, env::internal::Particle const& p2, vec3 const& r) const noexcept {
            const double r2 = r.norm_squared();
            if (cutoff_radius > 0.0 && r2 > cutoff_radius*cutoff_radius) return {};

            const double inv_r = 1.0 / std::sqrt(r2);
            double inv_pow = inv_r;
            for (int i = 0; i < exp; i++) inv_pow *= inv_r;
            const double mag = pre_factor * p1.mass * p2.mass * inv_pow;

            return mag * r;  // Force vector pointing along +r
        }

        [[nodiscard]] PowerLaw mix(PowerLaw const& other) const noexcept {
            // Arithmetic average of prefactor and cutoff
            const double mixed_prefactor = 0.5 * (pre_factor + other.pre_factor);
            const double mixed_cutoff = 0.5 * (cutoff_radius + other.cutoff_radius);
            return PowerLaw(static_cast<uint8_t>(std::round(mixed_prefactor)), mixed_cutoff);
        }
    };
}