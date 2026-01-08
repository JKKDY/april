#pragma once

#include <cstdint>
#include "april/math/vec3.hpp"


namespace april {
	// concepts
	template <typename T, typename... Ts>
	concept same_as_any = (... or std::same_as<T, Ts>);

	// vector aliases
	using vec3 = utils::Vec3<double>;
	using int3 = utils::Vec3<int32_t>;
	using uint3 = utils::Vec3<uint32_t>;

	// helpers
	using vec3_ptr = utils::Vec3Ptr<vec3::type>;

} // namespace april

