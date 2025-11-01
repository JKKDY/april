#pragma once

#include <string>
#include <type_traits>
#include <cstdint>
#include <sstream>
#include <any>

#include "april/common.h"


namespace april::env {
	enum class ParticleState : uint8_t {
		ALIVE      = 1u << 0, // Moves, exerts and experiences forces
		DEAD       = 1u << 1, // Inactive; no movement or interaction
		PASSIVE    = 1u << 2, // Moves, experiences forces but exerts none
		STATIONARY = 1u << 3, // Exerts forces but does not move or respond
		EXERTING   = ALIVE | STATIONARY, // Can exert forces on others
		MOVABLE    = ALIVE | PASSIVE,    // Can move (may or may not exert forces)
		ALL        = ~0u  // Matches all states
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


	using FieldMask = uint32_t;

	enum class Field : FieldMask {
		position     = 1u << 0,
		velocity     = 1u << 1,
		force        = 1u << 2,
		old_position = 1u << 3,
		old_force    = 1u << 4,
		state        = 1u << 5,
		mass         = 1u << 6,
		type         = 1u << 7,
		id           = 1u << 8,
		user_data    = 1u << 9,
		all			 = ~0u
	};

	constexpr FieldMask to_field_mask(Field f) { return static_cast<FieldMask>(f); }

	constexpr FieldMask operator|(const Field a, const Field b) { return to_field_mask(a) | to_field_mask(b); }
	constexpr FieldMask operator|(const FieldMask m, const Field f) { return m | to_field_mask(f); }
	constexpr FieldMask operator|(const Field a, const FieldMask m) { return to_field_mask(a) | m; }

	template<FieldMask M, Field F>
	inline constexpr bool has_field = (M & to_field_mask(F)) != 0;



	template< typename T, Field F, FieldMask M, bool Ref, bool Const = false>
	using make_field = std::conditional_t<Ref,
		std::conditional_t<Const, const T&, T&>,
		std::conditional_t<Const, const T, T>>;


	using ParticleType = uint16_t;
	using ParticleID = uint32_t;

	template <typename T>
	concept IsUserData =
		std::default_initializable<T> &&
		std::is_trivially_copyable_v<T> &&
		std::is_trivially_destructible_v<T> &&
		std::is_standard_layout_v<T> &&
		(!std::is_polymorphic_v<T>);

	struct NoUserData {};



	template<typename Data = NoUserData>
	struct ParticleData {
		using user_data_t = Data;
	};

	// user facing declaration with optional fields and non typed field for user data
	struct Particle {
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

		std::any user_data {}; // custom user data
	};


	// reference to particle data. This is passed to e.g. controllers & boundaries that can mutate particle data
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef {
		// everything by reference except for type & id
		make_field<vec3, Field::position, M, true> position;
		make_field<vec3, Field::velocity, M, true> velocity;
		make_field<vec3, Field::force,    M, true> force;

		make_field<vec3, Field::old_position, M, true> old_position;
		make_field<vec3, Field::old_force, M, true> old_force;

		make_field<double, Field::mass, M, true> mass;
		make_field<ParticleState, Field::state, M, true> state;
		make_field<ParticleType, Field::type, M, false, true> type; // const
		make_field<ParticleID, Field::id, M, false, true> id; // const

		make_field<UserDataT, Field::user_data, M, true> user_data;
	};

	// a restricted reference that only allows force to be mutated e.g. for fields
	template<FieldMask M, IsUserData UserDataT>
	struct RestrictedParticleRef {
		// everything except for numerical types by reference and const except for force
		make_field<vec3, Field::position, M, true, true> position;
		make_field<vec3, Field::velocity, M, true, true> velocity;
		make_field<vec3, Field::force,    M, true, false> force;

		make_field<vec3, Field::old_position, M, true, true> old_position;
		make_field<vec3, Field::old_force, M, true, true> old_force;

		make_field<double, Field::mass, M, false, true> mass;
		make_field<ParticleState, Field::state, M, false, true> state;
		make_field<ParticleType, Field::type, M, false, true> type;
		make_field<ParticleID, Field::id, M, false, true> id;

		make_field<UserDataT, Field::user_data, M, true, true> user_data;
	};


	// an immutable reference to the particle data. meant for read only e.g. in monitors
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView {
		// everything const
		make_field<vec3, Field::position, M, true, true> position;
		make_field<vec3, Field::velocity, M, true, true> velocity;
		make_field<vec3, Field::force,    M, true, true> force;

		make_field<vec3, Field::old_position, M, true, true> old_position;
		make_field<vec3, Field::old_force, M, true, true> old_force;

		make_field<double, Field::mass, M, false, true> mass;
		make_field<ParticleState, Field::state, M, false, true> state;
		make_field<ParticleType, Field::type, M, false, true> type;
		make_field<ParticleID, Field::id, M, false, true> id;

		make_field<UserDataT, Field::user_data, M, true, true> user_data;
	};


	// easy terminal diagnostics
	template<typename P>
	std::string particle_to_string(const P & p) {
		std::ostringstream oss;
		oss << "Particle ID: " << p.id << "\n"
			<< "Position: " << p.position.to_string() << "\n"
			<< "Velocity: " << p.velocity.to_string() << "\n"
			<< "Force: " << p.force.to_string() << "\n"
			<< "Mass: " << p.mass << "\n"
			<< "Type: " << p.type << "\n"
			<< "State: " << static_cast<int>(p.state) << "\n";
		return oss.str();
	}


	namespace internal
	{
		// used internally in system. Holds all data of a particle
		template<IsUserData ParticleData>
		struct ParticleRecord {

			ParticleRecord() = default;

			ParticleID id {};		// id of the particle.
			ParticleType type {};	// type of the particle.

			vec3 position;			// current position of the particle.
			vec3 old_position;		// previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;			// current velocity of the particle.
			vec3 force;				// current force acting on the particle.
			vec3 old_force;			// previous force acting on the particle.

			ParticleState state {};	// state of the particle.
			double mass {};			// mass of the particle.

			ParticleData user_data;

			bool operator==(const ParticleRecord& other) const {
				return id == other.id;
			}
		};
	}

}


