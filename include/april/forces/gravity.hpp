#pragma once

#include <cmath>

#include "april/common.hpp"
#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"


namespace april::force {
    struct Gravity : Force{
        static constexpr env::FieldMask fields = +env::Field::mass;

        double grav_constant;

        explicit Gravity(const double grav_const = 1.0, const double cutoff = no_cutoff)
            : Force(cutoff), grav_constant(grav_const) {}

        template<env::FieldMask M, env::IsUserData U>
        vec3 eval(const env::ParticleView<M, U> & p1, const env::ParticleView<M, U> & p2, const vec3& r) const noexcept {
            const double inv_r = 1.0 / r.norm();
            const double mag = grav_constant * p1.mass * p2.mass * inv_r * inv_r;

            return mag * inv_r * r;  // Force vector pointing along +r
        }

        [[nodiscard]] Gravity mix(Gravity const& other) const noexcept {
            // Arithmetic average of pre-factor and cutoff
            const double mixed_factor = 0.5 * (grav_constant + other.grav_constant);
            const double mixed_cutoff = 0.5 * (cutoff() + other.cutoff());
            return Gravity(static_cast<uint8_t>(std::round(mixed_factor)), mixed_cutoff);
        }

        bool operator==(const Gravity&) const = default;
    };
}