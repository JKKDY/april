#pragma once

#include <cstdint>
#include <april/base/bitmask.hpp>

namespace april::env {


	enum class Field : uint16_t {
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

	AP_ENABLE_BITMASK_OPERATORS(Field)


	constexpr Field to_field_mask(Field f) { return static_cast<Field>(f); }

	template<class T>
	concept HasFields = requires { std::remove_cvref_t<T>::fields; };

	template<HasFields Self>
	inline constexpr Field FieldOf = std::remove_cvref_t<Self>::fields;

	template<Field M, Field F>
	inline constexpr bool has_field_v = (M & F) != Field::none;
}




