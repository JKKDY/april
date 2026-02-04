#pragma once

#include "april/base/traits.hpp"

#include "april/env/data.hpp"
#include "april/particle/fields.hpp"
#include "april/particle/defs.hpp"

#include "april/forces/force.hpp"
#include "april/forces/force_table.hpp"

#include "april/boundaries/boundary.hpp"
#include "april/boundaries/boundary_table.hpp"

#include "april/controllers/controller.hpp"
#include "april/fields/field.hpp"
#include "april/utility/pack_storage.hpp"

namespace april::env::internal {
	template<class FPack, class BPack, class CPack, class FFPack, class U>
	   struct EnvironmentTraits;

	// this class holds relevant types derived from template parameter packs
	// this significantly cleans up dependent type declarations in other classes
	template<
	force::IsForce... Fs,
	boundary::IsBoundary... BCs,
	controller::IsController... Cs,
	field::IsField... FFs,
	IsUserData U
	>
	struct EnvironmentTraits<
		force::ForcePack<Fs...>,
		boundary::BoundaryPack<BCs...>,
		controller::ControllerPack<Cs...>,
		field::FieldPack<FFs...>,
		U
	>
	{
		// Core Packs
		using FPackT  = force::ForcePack<Fs...>;
		using BPack_t  = boundary::BoundaryPack<BCs...>;
		using CPack_t  = controller::ControllerPack<Cs...>;
		using FFPack_t = field::FieldPack<FFs...>;

		// Derived Variants
		using force_variant_t    = force::internal::VariantType_t<Fs...>;
		using boundary_variant_t = boundary::internal::VariantType_t<BCs...>;

		// Derived Storage Types
		using controller_storage_t = shared::internal::PackStorage<Cs...>;
		using field_storage_t      = shared::internal::PackStorage<FFs...>;

		// Table Types
		using boundary_table_t = boundary::internal::BoundaryTable<boundary_variant_t>;
		using force_table_t    = force::internal::ForceTable<force_variant_t>;

		// particles
		using user_data_t = U;
		using particle_record_t = ParticleRecord<user_data_t>;
		template<FieldMask M> using particle_ref_t = ParticleRef<M, user_data_t>;
		template<FieldMask M> using restricted_particle_ref_t = RestrictedParticleRef<M, user_data_t>;
		template<FieldMask M> using particle_view_t = ParticleView<M, user_data_t>;

		// Environment Data type
		using environment_data_t = EnvironmentData<
			force_variant_t,
			boundary_variant_t,
			controller_storage_t,
			field_storage_t>;

		// Validity check: check for membership in parameter packs
		template<typename T> static constexpr bool is_valid_force_v = same_as_any<T, Fs...>;
		template<typename T> static constexpr bool is_valid_boundary_v = same_as_any<T, BCs...>;
		template<typename T> static constexpr bool is_valid_controller_v = same_as_any<T, Cs...>;
		template<typename T> static constexpr bool is_valid_field_v = same_as_any<T, FFs...>;
	};


	// environment traits concept
	template<typename T>
    struct is_environment_traits : std::false_type {};

	template<typename... Packs>
	struct is_environment_traits<EnvironmentTraits<Packs...>> : std::true_type {};

	template<typename T>
	concept IsEnvironmentTraits = is_environment_traits<T>::value;



	// pack utilities for detecting and localizing parameter packs
	// base case
	template<template<class ... > class Pack, class ... Packs>
	struct contains_pack : std::false_type {};

	// recursion: skip until match found
	template<template<class ...> class Pack, class Head, class... Tail>
	struct contains_pack<Pack, Head, Tail...> : contains_pack<Pack, Tail...> {};

	// match: Head = Pack
	template<template<class ...> class Pack, class... Ts, class ... Tail>
	struct contains_pack<Pack, Pack<Ts...>, Tail...> : std::true_type {};

	// convenience alias
	template<template<class ...> class Pack, class ... Packs>
	static constexpr bool contains_pack_v = contains_pack<Pack, Packs...>::value;


	// base case: nothing found, return empty pack
	template<template<class ...> class Pack, class ... Packs>
	struct get_pack {using type = Pack<>; };

	// recursion: skip until non matches (Head != Pack<...>)
	template<template<class ...> class Pack, class Head, class ... Tail>
	struct get_pack<Pack, Head, Tail...> : get_pack<Pack, Tail...> {};

	// match pack: Head = Pack<...>
	template<template<class ...> class Pack, class... Ts, class... Tail>
	struct get_pack<Pack, Pack<Ts...>, Tail...> {
		static_assert(!contains_pack_v<Pack, Tail...>, "Duplicate pack of the same kind was provided.");
		using type = Pack<Ts...>;
	};

	template<template<class ...> class Pack, class ... Packs>
	using get_pack_t = get_pack<Pack, std::remove_cvref_t<Packs>...>::type;


	// TMP utilities for finding and locating user data type
	template<class T>
	struct is_particle_data : std::false_type {};

	template<class U>
	struct is_particle_data<ParticleData<U>> : std::true_type {};

	template<class T>
	inline constexpr bool is_particle_data_v = is_particle_data<std::remove_cvref_t<T>>::value;


	// Base Case
	template<class ... Args>
	struct get_user_data {using type = ParticleData<>::user_data_t;};

	// recursion
	template<class Head, class ... Tail>
	struct get_user_data<Head, Tail...> : get_user_data<Tail...>{};

	// match case:
	template<class UserData, class ... Tail>
	struct get_user_data<ParticleData<UserData>, Tail...> {
		static_assert(!(... || is_particle_data_v<Tail>), "Multiple Particle<...> markers provided to Environment().");
		using type = UserData;
	};

	template<class... Args>
	using get_user_data_t = get_user_data<std::remove_cvref_t<Args>...>::type;


	template<class T> static constexpr bool is_any_pack_v =
		force::IsForcePack<T> || boundary::IsBoundaryPack<T> ||
		controller::IsControllerPack<T> || field::IsFieldPack<T> ||
		is_particle_data_v<T>;

}
