#pragma once

#include <cstdint>

#include "april/shared/vec3.hpp"

/// Cross-compiler force-inline
#if defined(_MSC_VER)
#   define AP_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#   define AP_FORCE_INLINE inline __attribute__((always_inline))
#else
#   define AP_FORCE_INLINE inline
#endif


#if defined(_MSC_VER)
	// MSVC requires the vendor-specific attribute to guarantee layout optimization
	// because the standard attribute is often ignored for ABI stability.
	#define AP_NO_UNIQUE_ADDRESS [[msvc::no_unique_address]]
#else
	// GCC, Clang, and standard conformant compilers
	#define AP_NO_UNIQUE_ADDRESS [[no_unique_address]]
#endif


namespace april {
	using vec3 = utils::Vec3<double>;
	using int3 = utils::Vec3<int32_t>;
    using uint3 = utils::Vec3<uint32_t>;

	using vec3_ptr = utils::Vec3Ptr<vec3::type>;

	template <typename T, typename... Ts> concept same_as_any = (... or std::same_as<T, Ts>);


} // namespace april

