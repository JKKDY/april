#pragma once

#include <cstdint>
#include "april/math/vec3.hpp"


#ifndef VEC3_TYPE
#define VEC3_TYPE double
#endif


namespace april {
	// concepts
	template <typename T, typename... Ts>
	concept same_as_any = (... or std::same_as<T, Ts>);

	// vector aliases
	using vec3 = utils::Vec3<VEC3_TYPE>; // general purpose vec3 (e.g. particle data)
	using vec3f = utils::Vec3<float>; // float vec3
	using vec3d = utils::Vec3<double>; // double vec3


	using int3 = utils::Vec3<int32_t>;
	using uint3 = utils::Vec3<uint32_t>;

	// helpers
	using vec3_ptr = utils::Vec3Ptr<vec3::type>;

} // namespace april

