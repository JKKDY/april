#pragma once
#include "april/base/types.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"


namespace april::particle
{
    /**
     * @brief Materialized particle record passed to container implementations.
     *
     * ParticleRecord is part of APRIL's container extension API. It represents the
     * normalized particle data produced during build_system(...), before a container
     * stores the data in its chosen memory layout.
     */
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

        bool operator==(const ParticleRecord&) const = default;
    };
}
