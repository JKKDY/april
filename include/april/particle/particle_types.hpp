#pragma once

#include <cstdint>
#include <concepts>
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

	template<ParticleField F>
	constexpr const char* field_name = [] {
		if constexpr (F == ParticleField::position)     return "POSITION";
		if constexpr (F == ParticleField::velocity)     return "VELOCITY";
		if constexpr (F == ParticleField::force)        return "FORCE";
		if constexpr (F == ParticleField::old_position) return "OLD_POSITION";
		if constexpr (F == ParticleField::state)        return "STATE";
		if constexpr (F == ParticleField::mass)         return "MASS";
		if constexpr (F == ParticleField::type)         return "TYPE";
		if constexpr (F == ParticleField::id)           return "ID";
		if constexpr (F == ParticleField::attributes)    return "USER_DATA";
		return "UNKNOWN_FIELD";
	}();


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

	struct NoParticleAttributes {};

	namespace particle {
		template <typename T>
		concept IsParticleAttributes =
			std::default_initializable<T> &&
			std::is_trivially_copyable_v<T> &&
			std::is_trivially_destructible_v<T> &&
			std::is_standard_layout_v<T> &&
			(!std::is_polymorphic_v<T>);

		// template struct used to tell the environment what user data will be used
		template<typename Data = NoParticleAttributes>
		struct ParticleAttributes { using particle_attributes_t = Data; };
	}

	template<typename Data = NoParticleAttributes>
	inline constexpr particle::ParticleAttributes<Data> particle_attributes {};


	namespace particle::internal {
		template<class T>
		concept HasFields = requires { std::remove_cvref_t<T>::fields; };

		template<HasFields Self>
		inline constexpr ParticleField FieldOf = std::remove_cvref_t<Self>::fields;

		template<ParticleField M, ParticleField F>
		inline constexpr bool has_field_v = (M & F) != ParticleField::none;
	}
}












