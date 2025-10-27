#pragma once

#include <cmath>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::force {
    // Power force law force (e.g., gravity, Coulomb)

    struct PowerLaw : Force{
        double pre_factor; // e.g. G or Coulomb constant
        uint8_t exp;

        explicit PowerLaw(const uint8_t exp, const double pre_factor_ = 1.0, const double cutoff = -1.0)
            : Force(cutoff), pre_factor(pre_factor_), exp(exp) {}

        [[nodiscard]] vec3 eval(env::internal::Particle const& p1, env::internal::Particle const& p2, vec3 const& r) const noexcept {
            const double r2 = r.norm_squared();
            if (has_cutoff() && r2 > cutoff*cutoff) return {};

            const double inv_r = 1.0 / std::sqrt(r2);
            double inv_pow = inv_r;
            for (int i = 0; i < exp; i++) inv_pow *= inv_r;
            const double mag = pre_factor * p1.mass * p2.mass * inv_pow;

            return mag * r;  // Force vector pointing along +r
        }

        [[nodiscard]] PowerLaw mix(PowerLaw const& other) const noexcept {
            // Arithmetic average of pre-factor and cutoff
            const double mixed_factor = 0.5 * (pre_factor + other.pre_factor);
            const double mixed_cutoff = 0.5 * (cutoff + other.cutoff);
            return PowerLaw(static_cast<uint8_t>(std::round(mixed_factor)), mixed_cutoff);
        }
    };
}