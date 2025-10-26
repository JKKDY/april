#pragma once

#include <limits>
#include <string>
#include <type_traits>
#include <cstdint>

#include "april/common.h"

// TODO enforce proper particle interface so containers can store particles however they want
// TODO find a way to extend particles (maybe a void pointer)

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



	using ParticleType = int;
	using ParticleID = int;
	static_assert(std::is_same_v<ParticleID, ParticleType>);
	
	using ParticleTypePair = std::pair<ParticleType, ParticleType>;
	using ParticleIDPair = std::pair<ParticleID, ParticleID>;

	constexpr ParticleID PARTICLE_ID_DONT_CARE = std::numeric_limits<ParticleID>::min();

	struct Particle {
		// Particle() = default;
		// Particle(const vec3& position, const vec3& velocity, double mass, ParticleType type = 0, ParticleID id = PARTICLE_ID_UNDEFINED, ParticleState state = ParticleState::ALIVE);
		ParticleID id = PARTICLE_ID_DONT_CARE;  // The id of the particle.
		ParticleType type = 0;  				// The type of the particle.
		vec3 position;      					// The position of the particle.
		vec3 velocity;      					// The velocity of the particle.
		double mass{};        					// The mass of the particle.
		ParticleState state{};					// The state of the particle.
	};

	namespace internal
	{
		using ParticleType = unsigned int;
		using ParticleID = unsigned int;

		using ParticleTypePair = std::pair<ParticleType, ParticleType>;
		using ParticleIDPair = std::pair<ParticleID, ParticleID>;

		struct Particle {
			using State = ParticleState;

			Particle() = default;
			Particle(ParticleID id, const vec3& position, const vec3& velocity, double mass, ParticleType type,
				State state = State::ALIVE, const vec3& force = {}, const vec3& old_force = {}, const vec3& old_position = {});

			void update_position(const vec3& dx) noexcept;
			void update_velocity(const vec3& dv) noexcept;
			void update_force(const vec3& df) noexcept;
			void reset_force() noexcept;

			bool operator==(const Particle& other) const;

			vec3 position;			// current position of the particle.
			vec3 old_position;		// previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;			// current velocity of the particle.
			vec3 force;				// current force acting on the particle.
			vec3 old_force;			// previous force acting on the particle.

			State state;			// state of the particle.
			double mass{};			// mass of the particle.
			ParticleType type{};	// type of the particle.
			ParticleID id{};		// id of the particle.

			[[nodiscard]] std::string to_string() const;
		};
	}


	struct ParticleRef {
		explicit ParticleRef(internal::Particle & p);

		void update_position(const vec3& dx) const noexcept;
		void update_velocity(const vec3& dv) const noexcept;
		void update_force(const vec3& df) const noexcept;
		void reset_force() const noexcept;

		vec3 & position;
		vec3 & old_position;
		vec3 & velocity;
		vec3 & force;
		vec3 & old_force;

		ParticleState & state;

		double & mass;
		const internal::ParticleType & type;
		const internal::ParticleID & id;

		bool operator==(const internal::Particle& other) const;
		[[nodiscard]] std::string to_string() const;
	};

	struct RestrictedParticleRef {
		explicit RestrictedParticleRef(internal::Particle & p);

		const vec3 & position;
		const vec3 & velocity;
		vec3 & force;

		const ParticleState state;

		const double mass;
		const internal::ParticleType type;
		const internal::ParticleID id;

		[[nodiscard]] std::string to_string() const;
	};

	struct ParticleView {
		explicit ParticleView(const internal::Particle& p);

		const vec3& position;
		const vec3& old_position;
		const vec3& velocity;
		const vec3& force;
		const vec3& old_force;

		const internal::Particle::State&  state;
		const double&				  mass;
		const internal::ParticleType&     type;
		const internal::ParticleID&       id;

		bool operator==(const internal::Particle& other) const;
		[[nodiscard]] std::string to_string() const;
	};
}

