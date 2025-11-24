#pragma once
#include "april/common.h"
#include "april/particle/defs.h"


namespace april::env {

 // user facing declaration with optional fields and non typed field for user data
    struct Particle {
        std::optional<ParticleID> id;			// The id of the particle.
        ParticleType type = 0;  				// The type of the particle.

        vec3 position;      					// The position of the particle.
        vec3 velocity;      					// The velocity of the particle.

        double mass{};        					// The mass of the particle.
        ParticleState state = ParticleState::ALIVE;	// The state of the particle.

        // optional data e.g. if initializing from a simulation snapshot
        std::optional<vec3> old_position;		// previous position of the particle. Useful for applying boundary conditions
        std::optional<vec3> old_force;			// previous force acting on the particle.
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
        Particle& at(const double x, const double y, const double z) noexcept {
            position = {x,y,z}; return *this;
        }
        Particle& with_velocity(const vec3& v) noexcept {
            velocity = v; return *this;
        }
        Particle& with_velocity(const double x, const double y, const double z) noexcept {
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
        Particle& with_old_force(const vec3& v) noexcept {
            old_force = v; return *this;
        }
        Particle& with_force(const vec3& v) noexcept {
            force = v; return *this;
        }
        Particle& with_data(const std::any& v) noexcept {
            user_data = v; return *this;
        }
    };

}