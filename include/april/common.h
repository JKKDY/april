#pragma once

#include <cstdint>

#include "april/utils/vec3.hpp"



namespace april {
	using vec3 = utils::Vec3<double>;
	using int3 = utils::Vec3<int32_t>;
    using uint3 = utils::Vec3<uint32_t>;

	template <typename T, typename... Ts> concept same_as_any = (... or std::same_as<T, Ts>);


} // namespace april

