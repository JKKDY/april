#pragma once

#include <concepts>


namespace april::field {

	class Field {

	};

	// define controller concept
	template <class FFs>
	concept IsField = std::derived_from<FFs, Field>;

	// define controller Pack
	template<IsField...>
	struct FieldPack {};

	// constrained variable template
	template<class... FFs>
	requires (IsField<FFs> && ...)
	inline constexpr FieldPack<FFs...> fields {};
}