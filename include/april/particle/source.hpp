#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/fields.hpp"

namespace april::env {

	namespace internal {

		// A Poison struct
		// will throw at compile time if any operator call is instantiated
		template<ParticleField F>
	    struct AccessForbidden {
	        // Assignment & Compound Assignment
	        template<typename T> auto operator= (T&&) { trigger(); return *this; }
	        template<typename T> auto operator+=(T&&) { trigger(); return *this; }
	        template<typename T> auto operator-=(T&&) { trigger(); return *this; }
	        template<typename T> auto operator*=(T&&) { trigger(); return *this; }
	        template<typename T> auto operator/=(T&&) { trigger(); return *this; }
	        template<typename T> auto operator%=(T&&) { trigger(); return *this; }

	        // Increment & Decrement
	        auto operator++()    { trigger(); return *this; }
	        auto operator++(int) { trigger(); return *this; }
	        auto operator--()    { trigger(); return *this; }
	        auto operator--(int) { trigger(); return *this; }

	        // Binary Arithmetic
	        template<typename T> friend auto operator+(AccessForbidden, T&&) { trigger(); return AccessForbidden{}; }
	        template<typename T> friend auto operator-(AccessForbidden, T&&) { trigger(); return AccessForbidden{}; }
	        template<typename T> friend auto operator*(AccessForbidden, T&&) { trigger(); return AccessForbidden{}; }
	        template<typename T> friend auto operator/(AccessForbidden, T&&) { trigger(); return AccessForbidden{}; }

	        // Comparison
	        template<typename T> auto operator==(T&&) const { trigger(); return false; }
	        template<typename T> auto operator!=(T&&) const { trigger(); return false; }
	        template<typename T> auto operator< (T&&) const { trigger(); return false; }
	        template<typename T> auto operator> (T&&) const { trigger(); return false; }
	        template<typename T> auto operator<=(T&&) const { trigger(); return false; }
	        template<typename T> auto operator>=(T&&) const { trigger(); return false; }

	        // Logical Ops
	        auto operator!()  const { trigger(); return false; }
	        auto operator~()  const { trigger(); return *this; }
	        template<typename T> auto operator&(T&&)  const { trigger(); return *this; }
	        template<typename T> auto operator|(T&&)  const { trigger(); return *this; }
	        template<typename T> auto operator^(T&&)  const { trigger(); return *this; }

	        // Member Access & Dereference
	        auto& operator*()  const { trigger(); return *this; }
	        auto* operator->() const { trigger(); return this; }
	        template<typename T> auto operator[](T&&) const { trigger(); return *this; }

	        // Conversion (to catch logging/printing)
	        template<typename T> operator T() const { trigger(); return T{}; }

	        // Stream output (to catch std::cout << p.field)
	        friend std::ostream& operator<<(std::ostream& os, const AccessForbidden&) {
	            trigger();
	            return os;
	        }

	    private:
			// Note: in c++26 use std::format to print particle field in static assert
	        template<typename U = void>
	        static void trigger() {
	            static_assert(std::is_same_v<U, int>,
	                "\n\nError: Field Access Violation!\n"
	                "The requested field is not present in the current Accessor Mask.\n"
	                "Make sure that all used fields are in the particle mask\n");
	        }
	    };
	}


	// conditional switch between type and a void struct that will throw on access during compile time
	template<typename T, ParticleField F, ParticleField M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, internal::AccessForbidden<F>>;



	template<ParticleField M, IsUserData U, bool IsConst>
	struct ParticleSource {

		// selects T* or const T*
		template<typename T>
		using Ptr = std::conditional_t<IsConst, const T* AP_RESTRICT, T* AP_RESTRICT>;

		using Vec3PtrT = math::Vec3Ptr<std::conditional_t<IsConst, const vec3::type, vec3::type>>;

		// data pointers (optimized away if not in M)
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, ParticleField::force, M> force;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, ParticleField::position, M> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, ParticleField::velocity, M> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<Vec3PtrT, ParticleField::old_position, M> old_position;

		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<double>, ParticleField::mass, M> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleState>, ParticleField::state, M> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleType>, ParticleField::type, M> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<ParticleID>, ParticleField::id, M> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<Ptr<U>, ParticleField::user_data, M> user_data;

		// getter (Used by ParticleRef/View)
		template<ParticleField F>
		constexpr auto get() const noexcept {
			if constexpr (has_field_v<M, F>) {
				if constexpr (F == ParticleField::force) return force;
				else if constexpr (F == ParticleField::position) return position;
				else if constexpr (F == ParticleField::velocity) return velocity;
				else if constexpr (F == ParticleField::old_position) return old_position;
				else if constexpr (F == ParticleField::mass) return mass;
				else if constexpr (F == ParticleField::state) return state;
				else if constexpr (F == ParticleField::type) return type;
				else if constexpr (F == ParticleField::id) return id;
				else if constexpr (F == ParticleField::user_data) return user_data;
			} else {
				return static_cast<Ptr<void>>(nullptr);
			}
		}
	};
}



