#pragma once

#include <limits>
#include <string>
#include <type_traits>

#include "Common.h"

namespace april::env {
	class ParticleGrid;

	enum class ParticleState {
		ALIVE = 0x1,
		DEAD = 0x2,
		PASSIVE = 0x4, ///< will not exert forces on other particles
		STATIONARY = 0x8,
		ALL = 0xF
	};

	// // Overload bitwise OR (`|`)
	// inline unsigned int operator|(ParticleState a, ParticleState b) {
	// 	return static_cast<unsigned int>(a) | static_cast<unsigned int>(b);
	// }

	// // Overload bitwise AND (`&`)
	// inline ParticleState operator&(ParticleState a, ParticleState b) {
	// 	return static_cast<ParticleState>(
	// 		static_cast<unsigned int>(a) & static_cast<unsigned int>(b)
	// 	);
	// }

	// // Overload bitwise NOT (`~`)
	// inline ParticleState operator~(ParticleState a) {
	// 	return static_cast<ParticleState>(
	// 		~static_cast<unsigned int>(a)
	// 	);
	// }

	using ParticleType = int;
	using ParticleID = int;
	static_assert(std::is_same<ParticleID, ParticleType>::value);
	
	using ParticleTypePair = std::pair<ParticleType, ParticleType>;
	using ParticleIDPair = std::pair<ParticleID, ParticleID>;

	constexpr ParticleID PARTICLE_ID_UNDEFINED = std::numeric_limits<ParticleID>::min();

	struct Particle {
		ParticleID id = PARTICLE_ID_UNDEFINED;  // The id of the particle.
		ParticleType type = 0;  				// The type of the particle.
		vec3 position;      					// The position of the particle.
		vec3 velocity;      					// The velocity of the particle.
		double mass;        					// The mass of the particle.
		ParticleState state;					// The state of the particle.
	};

	namespace impl
	{
		using ParticleType = unsigned int;
		using ParticleID = unsigned int;

		using ParticleTypePair = std::pair<ParticleType, ParticleType>;
		using ParticleIDPair = std::pair<ParticleID, ParticleID>;

		struct Particle {
			using State = april::env::ParticleState;
			using ParticleIndex = size_t;

			Particle(size_t index, ParticleID id, const vec3& position, const vec3& velocity, double mass, ParticleType type,
				State state = State::ALIVE, const vec3& force = {}, const vec3& old_force = {}, const vec3& old_position = {});


			void update_position(const vec3& dx);
			void reset_force();

			bool operator==(const Particle& other) const;

			vec3 position;				// current position of the particle.
			vec3 old_position;			// previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;				// current velocity of the particle.
			vec3 force;					// current force acting on the particle.
			vec3 old_force;				// previous force acting on the particle.

			State state;				// state of the particle.

			const double mass;			// mass of the particle.
			const ParticleType type;    // type of the particle.
			const ParticleID id;		// id of the particle.
			const ParticleIndex index;	// index of the particle in the particle vector.

			std::string to_string() const;
		};
	}
}

