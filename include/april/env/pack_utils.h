#pragma once
#include <type_traits>

#include "april/forces/force.h"
#include "april/boundaries/boundary.h"
#include "april/controllers/controller.h"
#include "april/fields/field.h"


namespace april::env::internal {

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
	using get_pack_t = typename get_pack<Pack,  std::remove_cvref_t<Packs>...>::type;



	template<class T>
	static constexpr bool is_any_pack_v =
		force::IsForcePack<T> || boundary::IsBoundaryPack<T> ||
		controller::IsControllerPack<T> || field::IsFieldPack<T>;


}
