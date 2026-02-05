#pragma once
#include <cstdint>



#ifndef VEC3_TYPE
#define VEC3_TYPE double
#endif

#include "april/math/vec3.hpp"
#include "april/simd/packed.hpp"


namespace april {
	// vector aliases
	using vec3 = math::Vec3<VEC3_TYPE, double>; // general purpose vec3 (e.g. particle data)
	using vec3f = math::Vec3<float, double>; // float vec3
	using vec3d = math::Vec3<double>; // double vec3

	using packed = simd::Packed<VEC3_TYPE>;
	using packedd = simd::Packed<double>;
	using packedf = simd::Packed<float>;

	using pvec3 = math::Vec3<packed>;
	using pvec3f = math::Vec3<packedf>;
	using pvec3d = math::Vec3<packedf>;

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

