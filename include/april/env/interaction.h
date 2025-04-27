#pragma once
#include <concepts>
#include <functional>
#include <variant>
#include <math.h>


#include "april/env/particle.h"
#include "april/utils/map.hpp"
#include "common.h"


namespace april::env {

    // Base functor for force calculations
    struct Force {
        virtual ~Force() = default;

        // Main functor interface: computes force between p1 and p2
        virtual vec3 operator()(const impl::Particle& p1, const impl::Particle& p2, const vec3& r) const noexcept = 0;

        // Mixes this force with another force; returns a new heap-allocated Force
        virtual std::unique_ptr<Force> mix(const Force* other) const = 0;

        double cutoff_radius = -1;  // Negative value means no cutoff
    };


    struct NoForce final : Force {
        vec3 operator()(const impl::Particle&, const impl::Particle&, const vec3&) const noexcept override {
            return vec3{0.0, 0.0, 0.0};
        };

        std::unique_ptr<Force> mix(const Force* other) const override {
           return std::make_unique<NoForce>();
        }
    };


    // Lennard-Jones potential (12-6)
    struct LennardJones final : Force {
        LennardJones(const double epsilon, const double sigma, const double cutoff = -1)
            : epsilon(epsilon), sigma(sigma) {
            cutoff_radius = (cutoff < 0) ? 3.0 * sigma : cutoff;
        }

        vec3 operator()(const impl::Particle& p1, const impl::Particle& p2, const vec3& r) const noexcept override {
            const double r2 = r.norm_squared();
            if (cutoff_radius > 0 && r2 > cutoff_radius * cutoff_radius)
                return vec3{0.0, 0.0, 0.0};

            const double inv_r2 = 1.0 / r2;
            const double sigma_r2 = sigma * sigma * inv_r2;
            const double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
            const double sigma_r12 = sigma_r6 * sigma_r6;
            const double magnitude = 24.0 * epsilon * inv_r2 * (2.0 * sigma_r12 - sigma_r6);

            return magnitude * r;  // Force vector
        }

        std::unique_ptr<Force> mix(const Force* other) const override {
            const auto* o = dynamic_cast<const LennardJones*>(other);
            if (!o)
                throw std::invalid_argument("Cannot mix LennardJones with non-LennardJones force");
    
            // Lorentz-Berthelot mixing 
            double mixed_epsilon = std::sqrt(this->epsilon * o->epsilon);
            double mixed_sigma = 0.5* (this->sigma + o->sigma);
            double cutoff = std::sqrt(this->cutoff_radius * o->cutoff_radius);
            return std::make_unique<LennardJones>(mixed_epsilon, mixed_sigma, cutoff);
        }

        double epsilon;  // Depth of the potential well
        double sigma;    // Distance at which potential is zero
    };


    // Inverse-square law (e.g., gravity, Coulomb)
    struct InverseSquare final : Force {
        explicit InverseSquare(const double pre_factor = 1.0, const double cutoff = -1)
            : pre_factor(pre_factor) {
            cutoff_radius = cutoff;
        }

        vec3 operator()(const impl::Particle& p1, const impl::Particle& p2, const vec3& r) const noexcept override {
            const double distance = r.norm();
            if (cutoff_radius > 0 && distance > cutoff_radius)
                return vec3{0.0, 0.0, 0.0};

            const double magnitude = pre_factor * p1.mass * p2.mass / (distance * distance * distance);
            return -magnitude * r;  
        }

        std::unique_ptr<Force> mix(const Force* other) const override {
            const auto* o = dynamic_cast<const InverseSquare*>(other);
            if (!o)
                throw std::invalid_argument("Cannot mix InverseSquare with non-InverseSquare force");
            
            double k = 0.5 * (pre_factor + o->pre_factor); 
            double cutoff = 0.5 * (cutoff_radius + o->cutoff_radius);
            return std::make_unique<InverseSquare>(k, cutoff);
        }

        double pre_factor;  // G or k constant
    };


    // Harmonic spring force (Hooke's law)
    struct Harmonic final : Force {
        Harmonic(const double k, const double r0) : k(k), r0(r0) {}

        vec3 operator()(const impl::Particle& p1, const impl::Particle& p2, const vec3& r) const noexcept override {
            const double distance = r.norm();
            const double magnitude = k * (distance - r0) / distance;
            return -magnitude * r;  
        }

        std::unique_ptr<Force> mix(const Force* other) const override {
            return std::make_unique<NoForce>();
        }

        double k;   // Spring constant
        double r0;  // Equilibrium distance
    };


    template <typename T> concept IsForce = std::is_base_of_v<Force, T>;


    namespace impl {
        using ForcePtr = std::unique_ptr<Force>;


        template <typename T> concept ForceMap = requires(
            T t, size_t i, size_t j
        ) {
            {t.get(i,j) } -> std::same_as<Force*>;
            {t.key_size() } -> std::same_as<size_t>;
        };


        struct InteractionInfo {
            InteractionInfo(const bool pair_contains_types, std::pair<int, int> key_pair, ForcePtr force):
                pair_contains_types(pair_contains_types), force(std::move(force)) {
                if (key_pair.first < key_pair.second) {
                    this->key_pair = {key_pair.first, key_pair.second};
                } else {
                    this->key_pair = {key_pair.second, key_pair.first};
                }
            }
            bool pair_contains_types;
            std::pair<int, int> key_pair;
            ForcePtr force;
        };


        class InteractionManager {
            using TypeForceMap = utils::impl::DensePairMap<Force, ParticleType>;
            using IdForceMap = utils::impl::DensePairMap<Force, ParticleID>;

            static_assert(ForceMap<TypeForceMap>, "TypeForceMap must implement ForceMap interface");
            static_assert(ForceMap<IdForceMap>, "IdForceMap must implement ForceMap interface");

        public:
			InteractionManager() = default;
            void build(std::vector<InteractionInfo> & interaction_infos,
                const std::unordered_map<env::ParticleType, impl::ParticleType> & usr_types_to_impl_types,
                const std::unordered_map<env::ParticleID, impl::ParticleID> & usr_ids_to_impl_ids
            );

			vec3 evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const;
            
        private:
            TypeForceMap inter_type_forces;   // Forces between different particle types (e.g. type A ? type B)
            IdForceMap intra_particle_forces; // Forces between specific particle instances (by ID)
        };

	} // namespace impl

} // namespace april::env



