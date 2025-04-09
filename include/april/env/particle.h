#pragma once

#include <string>

#include "Common.h"

namespace april::env {
	class ParticleGrid;

	enum class ParticleState {
		ALIVE = 0x1,
		DEAD = 0x2,
		PASSIVE = 0x4, ///< will not exert forces on other particles
		STATIONARY = 0x8
	};

	struct Particle {
		Particle(int id, const vec3& position, const vec3& velocity, double mass, int type, ParticleState state) :
			id(id), position(position), velocity(velocity), mass(mass), type(type), state(state) {}

		int	id;             ///< The id of the particle.
		vec3 position;      ///< The position of the particle.
		vec3 velocity;      ///< The velocity of the particle.
		double mass;        ///< The mass of the particle.
		int type;           ///< The type of the particle.
		ParticleState state;///< The state of the particle.
	};

	namespace impl
	{

		struct Particle {
			using State = april::env::ParticleState;

			Particle(unsigned int id, const vec3& position, const vec3& velocity, double mass, unsigned int type,
				State state = State::ALIVE, const vec3& force = {}, const vec3& old_force = {}, const vec3& old_position = {});

			Particle(unsigned int id, const vec3& position, const vec3& velocity, double mass, unsigned int type, State state = State::ALIVE);


			void update_position(const vec3& dx);
			void reset_force();

			vec3 position;				///< The current position of the particle.
			vec3 old_position;			///< The previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;				///< The current velocity of the particle.
			vec3 force;					///< The current force acting on the particle.
			vec3 old_force;				///< The previous force acting on the particle.

			State state;				///< The state of the particle.

			const double mass;			///< The mass of the particle.
			const unsigned int type;    ///< The type of the particle.
			const unsigned int id;		///< The id of the particle.


			std::string to_string() const;
		};
	}
}

