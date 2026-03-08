#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/attributes.hpp"


namespace april::particle::internal {

	// A Poison struct
	// will throw at compile time if any operator call is instantiated
	template<ParticleField F>
    struct AccessForbidden {
        // Assignment & Compound Assignment
        template<typename T> auto operator= (T&&) const { trigger(); return *this; }
        template<typename T> auto operator+=(T&&) const { trigger(); return *this; }
        template<typename T> auto operator-=(T&&) const { trigger(); return *this; }
        template<typename T> auto operator*=(T&&) const { trigger(); return *this; }
        template<typename T> auto operator/=(T&&) const { trigger(); return *this; }
        template<typename T> auto operator%=(T&&) const { trigger(); return *this; }

        // Increment & Decrement
        auto operator++()    const { trigger(); return *this; }
        auto operator++(int) const { trigger(); return *this; }
        auto operator--()    const { trigger(); return *this; }
        auto operator--(int) const { trigger(); return *this; }

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
            static_assert(std::is_same_v<U, int>, // will evaluate to false
                "\n\nError: Field Access Violation!\n"
                "The requested field is not present in the current Accessor Mask.\n"
                "Make sure that all used fields are in the particle mask\n");
        }
    };


	// selects Mutable, Const, or Poison based on the masks and container constness
	template<typename MutableT, typename ConstT, ParticleField F, ParticleField ReadMask, ParticleField WriteMask>
	using field_access_t = std::conditional_t<
		has_field_v<ReadMask | WriteMask, F>, // is it accessible at all?
		std::conditional_t<
			has_field_v<WriteMask, F>,  // is it mutable?
			MutableT,
			ConstT
		>,
		AccessForbidden<F> // if not accessible return Poison
	>;


	template<ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
	struct ParticleSource {

		// type generator for standard scalar pointers
		template<typename T, ParticleField F>
		using Ptr = field_access_t<T* AP_RESTRICT, const T* AP_RESTRICT, F, ReadMask, WriteMask>;

		// type generator for Vec3 pointers
		template<ParticleField F>
		using Vec3PtrT = field_access_t<math::Vec3Ptr<vec3::type>, math::Vec3Ptr<const vec3::type>, F, ReadMask, WriteMask>;

		// data pointers (optimized away to empty poison structs if not accessible)
		AP_NO_UNIQUE_ADDRESS Vec3PtrT<ParticleField::force>        force;
		AP_NO_UNIQUE_ADDRESS Vec3PtrT<ParticleField::position>     position;
		AP_NO_UNIQUE_ADDRESS Vec3PtrT<ParticleField::velocity>     velocity;
		AP_NO_UNIQUE_ADDRESS Vec3PtrT<ParticleField::old_position> old_position;

		// scalar pointers (optimized away to empty poison structs if not accessible)
		AP_NO_UNIQUE_ADDRESS Ptr<double,        ParticleField::mass>       mass;
		AP_NO_UNIQUE_ADDRESS Ptr<ParticleState, ParticleField::state>      state;
		AP_NO_UNIQUE_ADDRESS Ptr<ParticleType,  ParticleField::type>       type;
		AP_NO_UNIQUE_ADDRESS Ptr<ParticleID,    ParticleField::id>         id;
		AP_NO_UNIQUE_ADDRESS Ptr<Attributes,    ParticleField::attributes> attributes;

		// ParticleField based getter
		template<ParticleField F>
		constexpr auto get() const noexcept {
			if constexpr (has_field_v<ReadMask | WriteMask, F>) {
				if constexpr (F == ParticleField::force) return force;
				else if constexpr (F == ParticleField::position) return position;
				else if constexpr (F == ParticleField::velocity) return velocity;
				else if constexpr (F == ParticleField::old_position) return old_position;
				else if constexpr (F == ParticleField::mass) return mass;
				else if constexpr (F == ParticleField::state) return state;
				else if constexpr (F == ParticleField::type) return type;
				else if constexpr (F == ParticleField::id) return id;
				else if constexpr (F == ParticleField::attributes) return attributes;
			} else {
				return AccessForbidden<F>{};
			}
		}
	};
}


