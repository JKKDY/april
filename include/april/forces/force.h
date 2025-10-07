#pragma once

#include <variant>
#include <cmath>
#include <concepts>

#include "april/common.h"
#include "april/env/particle.h"


namespace april::env {

    template<typename F> concept IsForce =
        std::copy_constructible<F> && std::assignable_from<F&, F const&> && std::movable<F> &&
    requires(F const& f,
             impl::Particle const& p1,
             impl::Particle const& p2,
             vec3 const& r,
             F const& o)
    {
        { f(p1, p2, r) } noexcept -> std::same_as<vec3>;
        { f.mix(o) } -> std::same_as<F>;
        { f.cutoff_radius } -> std::convertible_to<double>;

    };

    template<typename T>
    struct is_force_variant : std::false_type {};

    template<IsForce... Fs>
    struct is_force_variant<std::variant<Fs...>> : std::true_type {};

    template<typename T>
    concept ForceVariant = is_force_variant<T>::value;


    template<IsForce... Fs>
    struct ForcePack { /* empty */ };

    template<class... Fs>
    inline constexpr ForcePack<Fs...> forces{};


    // No-op force: always returns zero vector and mixes to itself.
    struct NoForce {
        // Negative cutoff_radius means "no cutoff"
        double cutoff_radius = 0.0;

        vec3 operator()(impl::Particle const&, impl::Particle const&, vec3 const&) const noexcept {
            return vec3{0.0, 0.0, 0.0};
        }

        [[nodiscard]] NoForce mix(NoForce const&) const noexcept {
            return {};
        }
    };


    // Lennard-Jones potential (12-6). epsilon: well depth; sigma: zero-cross distance.
    struct LennardJones {
        double epsilon; // Depth of the potential well
        double sigma; // Distance at which potential is zero
        double cutoff_radius; // Maximum interaction distance; negative: no cutoff

        LennardJones(const double epsilon_, const double sigma_, const double cutoff = -1.0)
        : epsilon(epsilon_), sigma(sigma_), sigma2(sigma * sigma) {
            cutoff_radius = (cutoff < 0.0) ? 3.0 * sigma : cutoff;
        }

        vec3 operator()(impl::Particle const&, impl::Particle const&, vec3 const& r) const noexcept {
            const double r2 = r.norm_squared();
            if (cutoff_radius > 0.0 && r2 > cutoff_radius * cutoff_radius)
                return vec3{0.0, 0.0, 0.0};

            const double inv_r2 = 1.0 / r2;
            const double sigma_r2 = sigma2 * inv_r2;
            const double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
            const double sigma_r12 = sigma_r6 * sigma_r6;
            const double magnitude = 24.0 * epsilon * inv_r2 * (2.0 * sigma_r12 - sigma_r6);

            // Force vector pointing along -r
            return -magnitude * r;
        }

        [[nodiscard]] LennardJones mix(LennardJones const& other) const noexcept {
            // Lorentz-Berthelot mixing rules
            const double mixed_epsilon = std::sqrt(epsilon * other.epsilon);
            const double mixed_sigma = 0.5 * (sigma + other.sigma);
            const double mixed_cutoff = std::sqrt(cutoff_radius * other.cutoff_radius);
            return {mixed_epsilon, mixed_sigma, mixed_cutoff};
        }
    private:
        double sigma2;
    };


    // Inverse-square law force (e.g., gravity, Coulomb). pre_factor scales G or k.
    struct InverseSquare {
        double pre_factor; // e.g. G or Coulomb constant
        double cutoff_radius; // Maximum interaction distance; negative: no cutoff

        explicit InverseSquare(const double pre_factor_ = 1.0, const double cutoff = -1.0)
            : pre_factor(pre_factor_), cutoff_radius(cutoff) {}

        vec3 operator()(impl::Particle const& p1, impl::Particle const& p2, vec3 const& r) const noexcept {
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

    // Harmonic spring force (Hooke's law). k: spring constant; r0: equilibrium distance.
    struct Harmonic {
        double k; // Spring constant
        double r0; // Equilibrium distance
        double cutoff_radius; // Negative: no cutoff

        Harmonic(const double k_, const double r0_)
        : k(k_), r0(r0_), cutoff_radius(-1.0) {}

        vec3 operator()(impl::Particle const&, impl::Particle const&, vec3 const& r) const noexcept {
            const double dist = r.norm();
            const double magnitude = k * (dist - r0) / dist; // F = -k * (dist - r0) * (r / dist)
            return -magnitude * r;
        }

        [[nodiscard]] Harmonic mix(Harmonic const& other) const noexcept {
            // Arithmetic average of k and r0; carry max cutoff
            const double mixed_k = 0.5 * (k + other.k);
            const double mixed_r0 = 0.5 * (r0 + other.r0);
            Harmonic h(mixed_k, mixed_r0);
            h.cutoff_radius = std::max(cutoff_radius, other.cutoff_radius);
            return h;
        }
    };


    namespace impl {

        template<ForceVariant FV> struct InteractionInfo {
            bool pair_contains_types;
            std::pair<int,int> key_pair;
            FV force;

            InteractionInfo(const bool is_type_pair, const std::pair<int,int>& key, FV f)
              : pair_contains_types(is_type_pair), key_pair(key), force(std::move(f))
            {}
        };

    }

} // namespace april::env