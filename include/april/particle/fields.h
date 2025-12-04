#pragma once

#include <type_traits>
#include <cstdint>


namespace april::env {

	using FieldMask = uint32_t;

	enum class Field : FieldMask {
		none			= 0u,
		position     	= 1u << 0,
		velocity     	= 1u << 1,
		force        	= 1u << 2,
		old_position 	= 1u << 3,
		old_force    	= 1u << 4,
		state        	= 1u << 5,
		mass         	= 1u << 6,
		type         	= 1u << 7,
		id           	= 1u << 8,
		user_data    	= 1u << 9,
		all			 	= ~0u
	};

	constexpr FieldMask to_field_mask(Field f) { return static_cast<FieldMask>(f); }
	constexpr FieldMask operator+(Field f) { return static_cast<FieldMask>(f); }

	constexpr FieldMask operator|(const Field a, const Field b) { return +a | +b; }
	constexpr FieldMask operator|(const FieldMask m, const Field f) { return m | +f;}
	constexpr FieldMask operator|(const Field a, const FieldMask m) { return +a | m; }


	template<class T>
	concept HasFields = requires { std::remove_cvref_t<T>::fields; };

	template<HasFields Self>
	inline constexpr FieldMask FieldOf = std::remove_cvref_t<Self>::fields;

	template<FieldMask M, Field F>
	inline constexpr bool has_field_v = (M & to_field_mask(F)) != 0;
}


