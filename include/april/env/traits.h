#pragma once

#include "april/env/data.h"

// #include "april/forces/force.h"
#include "april/forces/force_table.h"

// #include "april/boundaries/boundary.h"
#include "april/boundaries/boundary_table.h"

#include "april/controllers/controller.h"
#include "april/fields/field.h"
#include "april/shared/pack_storage.h"

namespace april::env::internal {
	template<class FPack, class BPack, class CPack, class FFPack>
	   struct EnvironmentTraits;


	template<
	force::IsForce... Fs,
	boundary::IsBoundary... BCs,
	controller::IsController... Cs,
	field::IsField... FFs
	>
	struct EnvironmentTraits<
		force::ForcePack<Fs...>,
		boundary::BoundaryPack<BCs...>,
		controller::ControllerPack<Cs...>,
		field::FieldPack<FFs...>
	>
	{
		// --- Core Packs ---
		using FPack_t  = force::ForcePack<Fs...>;
		using BPack_t  = boundary::BoundaryPack<BCs...>;
		using CPack_t  = controller::ControllerPack<Cs...>;
		using FFPack_t = field::FieldPack<FFs...>;

		// --- Derived Variants ---
		using force_variant_t    = std::variant<Fs...>;
		using boundary_variant_t = boundary::internal::VariantType_t<BCs...>;

		// --- Derived Storage Types ---
		using controller_storage_t = shared::internal::PackStorage<Cs...>;
		using field_storage_t      = shared::internal::PackStorage<FFs...>;

		// --- Table Types ---
		using boundary_table_t = boundary::internal::BoundaryTable<boundary_variant_t>;
		using force_table_t    = force::internal::ForceTable<force_variant_t>;

		// --- Environment Data type ---
		using environment_data_t = EnvironmentData<
			force_variant_t,
			boundary_variant_t,
			controller_storage_t,
			field_storage_t>;

		// This helper will be true if T is one of the types in Cs...
		template<typename T>
		static constexpr bool is_valid_force_v = same_as_any<T, Fs...>;

		// This helper will be true if T is one of the types in BCs...
		template<typename T>
		static constexpr bool is_valid_boundary_v = same_as_any<T, BCs...>;

		// This helper will be true if T is one of the types in Cs...
		template<typename T>
		static constexpr bool is_valid_controller_v = same_as_any<T, Cs...>;

		// This helper will be true if T is one of the types in FFs...
		template<typename T>
		static constexpr bool is_valid_field_v = same_as_any<T, FFs...>;
	};

	template<typename T>
    struct is_environment_traits : std::false_type {};

	template<typename... Packs>
	struct is_environment_traits<EnvironmentTraits<Packs...>> : std::true_type {};

	template<typename T>
	concept IsEnvironmentTraits = is_environment_traits<T>::value;
}
