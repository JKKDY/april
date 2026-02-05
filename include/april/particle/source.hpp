#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/defs.hpp"
#include "april/particle/fields.hpp"

namespace april::env {
	// conditional switch between type and monostate
	template<typename T, Field F, FieldMask M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, std::monostate>;



	template<FieldMask M, IsUserData U, bool IsConst>
	struct ParticleSource {

		// selects T* or const T*
		template<typename T>
		using Ptr = std::conditional_t<IsConst, const T* AP_RESTRICT, T* AP_RESTRICT>;

		using Vec3PtrT = math::Vec3Ptr<std::conditional_t<IsConst, const vec3::type, vec3::type>>;

		// data pointers (optimized away if not in M)
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, Field::force, M> force;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, Field::position, M> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, Field::velocity, M> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, Field::old_position, M> old_position;

		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<double>, Field::mass, M> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleState>, Field::state, M> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleType>, Field::type, M> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleID>, Field::id, M> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<U>, Field::user_data, M> user_data;

		// getter (Used by ParticleRef/View)
		template<Field F>
		constexpr auto get() const noexcept {
			if constexpr (has_field_v<M, F>) {
				if constexpr (F == Field::force) return force;
				else if constexpr (F == Field::position) return position;
				else if constexpr (F == Field::velocity) return velocity;
				else if constexpr (F == Field::old_position) return old_position;
				else if constexpr (F == Field::mass) return mass;
				else if constexpr (F == Field::state) return state;
				else if constexpr (F == Field::type) return type;
				else if constexpr (F == Field::id) return id;
				else if constexpr (F == Field::user_data) return user_data;
			} else {
				return static_cast<Ptr<void>>(nullptr);
			}
		}
	};
}