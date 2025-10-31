#pragma once

#include <string>
#include <type_traits>
#include <cstdint>
#include <sstream>

#include "april/common.h"


namespace april::env {
	enum class ParticleState : uint8_t {
		ALIVE      = 0b00000001, // Moves, exerts and experiences forces
		DEAD       = 0b00000010, // Inactive; no movement or interaction
		PASSIVE    = 0b00000100, // Moves, experiences forces but exerts none
		STATIONARY = 0b00001000, // Exerts forces but does not move or respond
		EXERTING   = ALIVE | STATIONARY, // Can exert forces on others
		MOVABLE    = ALIVE | PASSIVE,    // Can move (may or may not exert forces)
		ALL        = 0b11111111  // Matches all states
	};


	// Bitwise OR
	inline ParticleState operator|(ParticleState a, ParticleState b) {
		return static_cast<ParticleState>(
			static_cast<unsigned int>(a) | static_cast<unsigned int>(b)
		);
	}

	// Bitwise AND
	inline ParticleState operator&(ParticleState a, ParticleState b) {
		return static_cast<ParticleState>(
			static_cast<unsigned int>(a) & static_cast<unsigned int>(b)
		);
	}

	// Bitwise NOT
	inline ParticleState operator~(ParticleState a) {
		return static_cast<ParticleState>(
			~static_cast<unsigned int>(a)
		);
	}

	// OR-assignment
	inline ParticleState& operator|=(ParticleState& a, ParticleState b) {
		a = a | b;
		return a;
	}

	// AND-assignment
	inline ParticleState& operator&=(ParticleState& a, ParticleState b) {
		a = a & b;
		return a;
	}


	template <typename T>
	concept IsUserData =
		std::default_initializable<T> &&
		std::is_trivially_copyable_v<T> &&
		std::is_trivially_destructible_v<T> &&
		std::is_standard_layout_v<T> &&
		(!std::is_polymorphic_v<T>);

	struct NoUserData {};

	using ParticleType = uint16_t;
	using ParticleID = uint32_t;


	template<IsUserData UserDataT = NoUserData>
	struct Particle {
		using user_data_t = UserDataT;

		std::optional<ParticleID> id;			// The id of the particle.
		ParticleType type = 0;  				// The type of the particle.

		vec3 position;      					// The position of the particle.
		vec3 velocity;      					// The velocity of the particle.

		double mass{};        					// The mass of the particle.
		ParticleState state{};					// The state of the particle.

		// optional data e.g. if initializing from a simulation snapshot
		std::optional<vec3> old_position = {};		// previous position of the particle. Useful for applying boundary conditions
		std::optional<vec3> old_force = {};			// previous force acting on the particle.
		std::optional<vec3> force = {};				// current force

		UserDataT user_data {}; // custom user data
	};


	template<IsUserData UserDataT>
	struct ParticleRef {

		vec3 & position;
		vec3 & velocity;
		vec3 & force;

		vec3 & old_position;
		vec3 & old_force;

		ParticleState & state;

		double & mass;
		const ParticleType type;
		const ParticleID id;

		UserDataT & user_data;

		bool operator==(const ParticleRef & other) const noexcept {return id == other.id; };
	};

	template<IsUserData UserDataT>
	struct RestrictedParticleRef {
		const vec3 & position;
		const vec3 & velocity;
		vec3 & force;

		const ParticleState state;

		const double mass;
		const ParticleType type;
		const ParticleID id;

		UserDataT & user_data;
	};

	template<IsUserData UserDataT>
	struct ParticleView {

		const vec3& position;
		const vec3& old_position;
		const vec3& velocity;
		const vec3& force;
		const vec3& old_force;

		const ParticleState  state;
		const double		 mass;
		const ParticleType   type;
		const ParticleID     id;

		const UserDataT & user_data;
	};


	template<typename P>
	std::string particle_to_string(const P & p) {
		std::ostringstream oss;
		oss << "Particle ID: " << p.id << "\n"
			<< "Position: " << p.position.to_string() << "\n"
			<< "Velocity: " << p.velocity.to_string() << "\n"
			<< "Force: " << p.force.to_string() << "\n"
			<< "Mass: " << p.mass << "\n"
			<< "Type: " << p.pe << "\n"
			<< "State: " << static_cast<int>(p.state) << "\n";
		return oss.str();
	}


	namespace internal
	{
		template<IsUserData UserDataT>
		struct Particle {

			Particle() = default;

			ParticleID id {};		// id of the particle.
			ParticleType type {};	// type of the particle.

			vec3 position;			// current position of the particle.
			vec3 old_position;		// previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;			// current velocity of the particle.
			vec3 force;				// current force acting on the particle.
			vec3 old_force;			// previous force acting on the particle.

			ParticleState state {};	// state of the particle.
			double mass {};			// mass of the particle.

			UserDataT user_data;

			bool operator==(const Particle& other) const {
				return id == other.id;
			}
		};
	}

}


