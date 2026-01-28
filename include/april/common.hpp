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
	using vec3 = math::Vec3<VEC3_TYPE>; // general purpose vec3 (e.g. particle data)
	using vec3f = math::Vec3<float>; // float vec3
	using vec3d = math::Vec3<double>; // double vec3


	using int3 = math::Vec3<int32_t>;
	using uint3 = math::Vec3<uint32_t>;

	// helpers
	using vec3_ptr = math::Vec3Ptr<vec3::type>;


	// enums
	enum class ExecutionPolicy : uint8_t {
		Seq = 0,        // Serial (Main thread only)
		Par = 1,        // Parallel (OMP Parallel For)
		Simd = 2,       // Vectorized (OMP Simd)
		ParSimd = 3     // Parallel + Vectorized
	};

} // namespace april

