/**
* @file source.hpp
 * @brief Compile-time field access control and zero-overhead data source abstraction.
 *
 * This file implements APRIL's "Poisoning" system.
 *
 * Any attempt to access a particle field that was not explicitly declared in a kernel's
 * Read/Write mask results in a clean compile-time error instead of a runtime bug.
 *
 * Forbidden fields are replaced with the AccessForbidden struct, which deliberately
 * triggers a static_assert on any operator usage. This gives excellent error messages
 * while maintaining zero runtime overhead.
 */
#pragma once

#include <ostream>
#include <type_traits>

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"


namespace april::particle::internal {

	/**
	 * Poison struct used to occupy slots for forbidden particle fields.
	 * Triggers a hard compiler error if any operator call is instantiated.
	 * This effectively "poisons" the field, preventing unauthorized usage.
	 */
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
            static_assert(std::is_same_v<U, int>, // will evaluate to false. One could also directly type false but clion/clangd complains
                "\n\nError: Field Access Violation!\n"
                "The requested field is not present in the current Accessor Mask.\n"
                "Make sure that all used fields are in the particle mask\n");
        }
    };


	/**
	 * Selects the correct type for a field based on the Read/Write masks:
	 *   - Mutable pointer/reference   if in WriteMask
	 *   - Const pointer/reference     if only in ReadMask
	 *   - AccessForbidden<F>          if not requested at all
	 */
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


	template<
		template <ParticleField> typename get_field_ptr
	>
	struct ParticleSourceStorageTypes {
		using force = std::invoke_result_t<get_field_ptr<ParticleField::force>, size_t>;
		using position = std::invoke_result_t<get_field_ptr<ParticleField::position>, size_t>;

	};


	/**
	 * Lightweight source of raw pointers to particle data.
	 *
	 * This struct acts as the common foundation for both ScalarParticleRef and
	 * PackedParticleRef. Thanks to APRIL_NO_UNIQUE_ADDRESS and the poisoning system,
	 * fields that are not requested compile down to zero bytes with no overhead.
	 *
	 * The masks are supplied by APRIL's kernel wrappers through their static
	 * `Read` and `Write` members.
	 * */
	template<
		ParticleField ReadMask,
		ParticleField WriteMask,
		typename Getter
	>
	struct ParticleSource {
	private:
		template<ParticleField F>
		using inferred_field_t = std::remove_cvref_t<
			decltype(std::declval<Getter&>().template operator()<F>())
		>;

		// lazy compile time conditional
		template<bool Enabled, ParticleField F>
		struct field_type {
			using type = AccessForbidden<F>;
		};

		template<ParticleField F>
		struct field_type<true, F> {
			using type = inferred_field_t<F>;
		};

		template<ParticleField F>
		using field_t = typename field_type<
			has_field_v<ReadMask | WriteMask, F>,
			F
		>::type;

		template<ParticleField F>
		static constexpr auto init_field(Getter& getter) {
			if constexpr (has_field_v<ReadMask | WriteMask, F>) {
				return getter.template operator()<F>();
			} else {
				return AccessForbidden<F>{};
			}
		}

	public:
		explicit ParticleSource(Getter& getter)
			: force(init_field<ParticleField::force>(getter))
			, position(init_field<ParticleField::position>(getter))
			, velocity(init_field<ParticleField::velocity>(getter))
			, old_position(init_field<ParticleField::old_position>(getter))
			, mass(init_field<ParticleField::mass>(getter))
			, state(init_field<ParticleField::state>(getter))
			, type(init_field<ParticleField::type>(getter))
			, id(init_field<ParticleField::id>(getter))
			, attributes(init_field<ParticleField::attributes>(getter))
		{}

		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::force> force;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::position> position;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::velocity> velocity;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::old_position> old_position;

		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::mass> mass;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::state> state;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::type> type;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::id> id;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::attributes> attributes;

		template<ParticleField F>
		constexpr decltype(auto) get() const noexcept {
			if constexpr (has_field_v<ReadMask | WriteMask, F>) {
				if constexpr (F == ParticleField::force) return (force);
				else if constexpr (F == ParticleField::position) return (position);
				else if constexpr (F == ParticleField::velocity) return (velocity);
				else if constexpr (F == ParticleField::old_position) return (old_position);
				else if constexpr (F == ParticleField::mass) return (mass);
				else if constexpr (F == ParticleField::state) return (state);
				else if constexpr (F == ParticleField::type) return (type);
				else if constexpr (F == ParticleField::id) return (id);
				else if constexpr (F == ParticleField::attributes) return (attributes);
			} else {
				return AccessForbidden<F>{};
			}
		}
	};

	template<ParticleField ReadMask, ParticleField WriteMask, typename Getter>
	auto make_particle_source(Getter& getter) {
		return ParticleSource<
			ReadMask,
			WriteMask,
			std::remove_cvref_t<Getter>
		>(getter);
	}
}
