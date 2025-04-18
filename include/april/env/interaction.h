#pragma once
#include <functional>
#include <variant>


#include <ankerl/unordered_dense.h>

#include "april/env/particle.h"
#include "common.h"


namespace april::env {

    // Base functor for force calculations
    struct Force {
        virtual ~Force() = default;

        // Main functor interface: computes force between p1 and p2
        virtual vec3 operator()(const Particle& p1, const Particle& p2, const vec3& r) const = 0;

        double cutoff_radius = -1;  // Negative value means no cutoff
    };

    struct NoForce : Force {
        virtual vec3 operator()(const Particle&, const Particle&, const vec3&) const override {
            return vec3{0.0, 0.0, 0.0};
        };
    };

    // Lennard-Jones potential (12-6)
    struct LennardJones : Force {
        LennardJones(double epsilon, double sigma, double cutoff = -1) 
            : epsilon(epsilon), sigma(sigma) {
            cutoff_radius = (cutoff < 0) ? 3.0 * sigma : cutoff;
        }

        vec3 operator()(const Particle& p1, const Particle& p2, const vec3& r) const override {
            const double r2 = r.norm_squared();
            if (cutoff_radius > 0 && r2 > cutoff_radius * cutoff_radius) {
                return vec3{0.0, 0.0, 0.0};
            }

            const double inv_r2 = 1.0 / r2;
            const double sigma_r2 = sigma * sigma * inv_r2;
            const double sigma_r6 = sigma_r2 * sigma_r2 * sigma_r2;
            const double sigma_r12 = sigma_r6 * sigma_r6;
            const double magnitude = 24.0 * epsilon * inv_r2 * (2.0 * sigma_r12 - sigma_r6);

            return magnitude * r;  // Force vector
        }

        double epsilon;  // Depth of the potential well
        double sigma;    // Distance at which potential is zero
    };

    // Inverse-square law (e.g., gravity, Coulomb)
    struct InverseSquare : Force {
        InverseSquare(double pre_factor = 1.0, double cutoff = -1) 
            : pre_factor(pre_factor) {
            cutoff_radius = cutoff;
        }

        vec3 operator()(const Particle& p1, const Particle& p2, const vec3& r) const override {
            const double distance = r.norm();
            if (cutoff_radius > 0 && distance > cutoff_radius) {
                return vec3{0.0, 0.0, 0.0};
            }

            const double magnitude = pre_factor * p1.mass * p2.mass / (distance * distance * distance);
            return -magnitude * r;  // Attractive force
        }

        double pre_factor;  // G or k constant
    };

    // Harmonic spring force (Hooke's law)
    struct Harmonic : Force {
        Harmonic(double k, double r0) : k(k), r0(r0) {}

        vec3 operator()(const Particle& p1, const Particle& p2, const vec3& r) const override {
            const double distance = r.norm();
            const double magnitude = k * (distance - r0) / distance;
            return -magnitude * r;  // Restoring force
        }

        double k;   // Spring constant
        double r0;  // Equilibrium distance
    };

    template <typename T> concept IsForce = std::is_base_of_v<Force, T>;

    struct Interaction {
        bool pair_contains_types;
        std::pair<int, int> id_or_type_pair;
        std::unique_ptr<Force> f;
    };

    namespace impl {


        class InteractionManager {
            
            using ForcePtr = std::unique_ptr<Force>;
			using Forces = std::vector<std::unique_ptr<Force>>;

        public:
			InteractionManager();
            void build(const std::vector<Interaction> & interactions);

			vec3 evaluate(const Particle& p1, const Particle& p2, const vec3& distance) const;
            
        private:
			std::unique_ptr<Force> mix_forces(ForcePtr force1, ForcePtr force2);

            std::vector<ForcePtr> inter_type_forces;   ///< Forces between different particle types (e.g. type A ? type B)
            std::vector<ForcePtr> intra_particle_forces; ///< Forces between specific particle instances (by ID)
        };

	} // namespace impl

} // namespace april::env



