#pragma once

#include <any>

#include "april/base/types.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/attributes.hpp"

namespace april {

    /**
     * @brief User-facing particle declaration used during environment setup.
     *
     * Particle is a declarative input object. During build_system(...), particles
     * are converted into the container-specific storage representation used by the
     * simulation System.
     */
    struct Particle {
        std::optional<ParticleID> id;        ///< Optional user-defined particle ID.
        ParticleType type = 0;               ///< User-defined particle type used for interaction lookup.

        vec3 position{};                     ///< Initial position.
        vec3 velocity{};                     ///< Initial velocity.

        double mass = 1.0;                   ///< Particle mass.
        ParticleState state = ParticleState::ALIVE;	// The state of the particle.

        // optional data e.g. if initializing from a simulation snapshot
        std::optional<vec3> old_position;		// previous position of the particle. Useful for applying boundary conditions
        std::optional<vec3> force;				// current force

        std::any user_data {}; // custom user data

        Particle& with_id(ParticleID v) noexcept {
            id = v; return *this;
        }
        Particle& as_type(const ParticleType v) noexcept {
            type = v; return *this;
        }
        Particle& at(const vec3& v) noexcept {
            position = v; return *this;
        }
        Particle& at(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            position = {x,y,z}; return *this;
        }
        Particle& with_velocity(const vec3& v) noexcept {
            velocity = v; return *this;
        }
        Particle& with_velocity(const vec3::type x, const vec3::type y, const vec3::type z) noexcept {
            velocity = {x,y,z}; return *this;
        }
        Particle& with_mass(const double v) noexcept {
            mass = v; return *this;
        }
        Particle& with_state(const ParticleState v) noexcept {
            state = v; return *this;
        }
        Particle& with_old_position(const vec3& v) noexcept {
            old_position = v; return *this;
        }
        Particle& with_force(const vec3& v) noexcept {
            force = v; return *this;
        }
        Particle& with_data(std::any v) {
            user_data = std::move(v);
            return *this;
        }
    };


    namespace particle
    {
        // used internally in system. Holds all data of a particle
        template<IsParticleAttributes A>
        struct ParticleRecord {
            using particle_attributes_t = A;
            ParticleRecord() = default;

            vec3 position;			// current position of the particle.
            vec3 force;				// current force acting on the particle.
            vec3 velocity;			// current velocity of the particle.
            vec3 old_position;		// previous position of the particle. Useful for applying boundary conditions

            double mass {};			// mass of the particle.
            ParticleState state {};	// state of the particle.
            ParticleID id {};		// id of the particle.
            ParticleType type {};	// type of the particle.

            [[no_unique_address]] A attributes; // optional user data

            bool operator==(const ParticleRecord& other) const {
                return id == other.id;
            }
        };
    }
}














