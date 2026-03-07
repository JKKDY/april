#pragma once

#include <cstdint>
#include <type_traits>

#include "april/base/bitmask.hpp"

namespace april {

	enum class ParticleField : uint16_t {
		none			= 0u,
		position     	= 1u << 0,
		velocity     	= 1u << 1,
		force        	= 1u << 2,
		old_position 	= 1u << 3,
		state        	= 1u << 4,
		mass         	= 1u << 5,
		type         	= 1u << 6,
		id           	= 1u << 7,
		attributes    	= 1u << 8,
		all			 	= static_cast<uint16_t>(~0u)
	};
	AP_ENABLE_BITMASK_OPERATORS(ParticleField)


	enum class ParticleState : uint8_t {
		ALIVE      = 1u << 0, // Moves, exerts and experiences forces
		DEAD       = 1u << 1, // Inactive; no movement or interaction
		PASSIVE    = 1u << 2, // Moves, experiences forces but exerts none
		STATIONARY = 1u << 3, // Exerts forces but does not move or respond
		INVALID	   = 1u << 7, // a special sentinel state to mark invalid data (e.g. to mark gaps in memory)

		EXERTING   = ALIVE | STATIONARY, // Can exert forces on others
		MOVABLE    = ALIVE | PASSIVE,    // Can move (may or may not exert forces)

		ALL        = 0b01111111  // Matches all states (Except invalid)
	};
	AP_ENABLE_BITMASK_OPERATORS(ParticleState)

	using ParticleType = uint16_t;
	using ParticleID = uint32_t;


	namespace particle::internal {
		template<class T>
		concept HasFields = requires { std::remove_cvref_t<T>::fields; };

		template<HasFields Self>
		inline constexpr ParticleField FieldOf = std::remove_cvref_t<Self>::fields;

		template<ParticleField M, ParticleField F>
		inline constexpr bool has_field_v = (M & F) != ParticleField::none;
	}
}
