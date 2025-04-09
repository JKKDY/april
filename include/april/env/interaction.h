#pragma once
#include <functional>

#include "april/env/particle.h"
#include "common.h"


namespace april::env {

    using ForceType = std::variant<LennardJones, InverseSquare, Harmonic>;

    struct ForceBase {
        double cutoff = 0;
    };

    /**
     * @brief Represents the Lennard-Jones force between two particles.
     */
    struct LennardJones : ForceBase {
        /**
         * @brief Default constructor.
         */
        LennardJones() = default;

        /**
         * @brief Constructs a Lennard-Jones object with parameters.
         * @param epsilon The lennard-jones epsilon.
         * @param sigma The lennard-jones sigma.
         * @param cutoff The cutoff-radius.
         */
        LennardJones(const double epsilon, const double sigma, const double cutoff)
            : ForceBase(sigma <= 0 ? cutoff = 3 * sigma), epsilon(epsilon), sigma(sigma) {
        }

        double epsilon;
        double sigma;
    };

    /**
     * @brief Represents an inverse-square force between two particles.
     */
    struct InverseSquare : ForceBase {
        /**
         * @brief Default constructor.
         */
        InverseSquare() = default;

        /**
         * @brief Constructs an inverse-square force object with parameters.
         * @param pre_factor The pre-factor of inverse square force.
         * @param cutoff The cutoff-radius.
         */
        InverseSquare(const double pre_factor=1, const double cutoff)
            : ForceBase(cutoff), pre_factor(pre_factor) {
        }

        double pre_factor;
    };


    struct Harmonic : ForceBase {
        Harmonic() = default;

        /**
        * @brief Harmonic force for membrane.
        */
        Harmonic(const double k, const double r0, const double cutoff)
            : ForceBase(cutoff), k(k), r0(r0) {
        }

        double k;
        double r0;
    };

    namespace impl {

        class Force {
        public:
            using ForceFunc = std::function<vec3(const vec3&, const Particle&, const Particle&)>;

            Force();

            explicit Force(ForceFunc force_func, double cutoff = NO_FORCE_CUTOFF);
            vec3 operator()(const vec3& diff, const Particle& p1, const Particle& p2) const;

            double cutoff() const;

        private:
            double cutoff_radius;   ///< The cutoff radius for the calculations.
            ForceFunc force_func{}; ///< The force function used for the calculations.
        };


        class InteractionManager {
            using ParticleType = unsigned int;
            using ParticleID = unsigned int;
            using ParticleTypePair = std::pair<ParticleType, ParticleType>;
            using ParticleIDPair = std::pair<ParticleID, ParticleID>;
            using GlobalForces = ankerl::unordered_dense::map<ParticleTypePair, Force, ForceKeyHash>;
            using LocalizedForces = ankerl::unordered_dense::map<ParticleIDPair, Force, ForceKeyHash>;


            struct ForceKeyHash {
                template <typename T1, typename T2>
                std::size_t operator()(const std::pair<T1, T2>& key) const {
                    return std::hash<T1>()(key.first) ^ (std::hash<T2>()(key.second) << 1);
                }
            };

        public:
            InteractionManager();

            void init();

            void add_force(const ForceType& force, int particle_type);
            void add_force(const ForceType& force, const ParticleIDPair& particle_ids);

            vec3 evaluate(const vec3& diff, const Particle& p1, const Particle& p2) const;

            double cutoff() const;

        private:

            static Force mix_forces(const ForceType& force1, const ForceType& force2);

            ankerl::unordered_dense::map<ParticleTypePair, Force, ForceKeyHash> global_forces; ///< forces between particle types.
            ankerl::unordered_dense::map<ParticleIDPair, Force, ForceKeyHash> localized_forces; ///< forces between specific particles.

            double cutoff_radius; ///< The cutoff radius.
        };
	} // namespace impl

} // namespace april::env



