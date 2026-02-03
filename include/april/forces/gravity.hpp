#pragma once

#include <cmath>

#include "april/base/types.hpp"
#include "april/particle/fields.hpp"
#include "april/forces/force.hpp"


namespace april::force {
    struct Gravity : Force{
        static constexpr env::FieldMask fields = +env::Field::mass;

        double grav_constant;

        explicit Gravity(const double grav_const = 1.0, const double cutoff = no_cutoff)
            : Force(cutoff), grav_constant(grav_const) {}

        template<typename vec3_t>
        vec3_t eval(const auto & p1, const auto & p2, const vec3_t& r) const noexcept {
            const double inv_r = 1.0 / r.norm();
            const double mag = grav_constant * p1.mass * p2.mass * inv_r * inv_r;

            return mag * inv_r * r;  // Force vector pointing along +r
        }

        [[nodiscard]] Gravity mix(Gravity const& other) const {
            if (std::abs(grav_constant - other.grav_constant) > 1e-9) {
                throw std::invalid_argument("Cannot mix different Gravitational Constants!");
            }
            return Gravity(grav_constant, std::max(cutoff(), other.cutoff()));
        }

        bool operator==(const Gravity&) const = default;
    };
}