#pragma once

#include <cstdint>
#include <april/base/bitmask.hpp>

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
		user_data    	= 1u << 8,
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
		if constexpr (F == ParticleField::user_data)    return "USER_DATA";
		return "UNKNOWN_FIELD";
	}();

	namespace env {
		template<class T>
		concept HasFields = requires { std::remove_cvref_t<T>::fields; };

		template<HasFields Self>
		inline constexpr ParticleField FieldOf = std::remove_cvref_t<Self>::fields;

		template<ParticleField M, ParticleField F>
		inline constexpr bool has_field_v = (M & F) != ParticleField::none;
	}
}






