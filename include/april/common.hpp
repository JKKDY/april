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


namespace april {
	using vec3 = utils::Vec3<double>;
	using int3 = utils::Vec3<int32_t>;
    using uint3 = utils::Vec3<uint32_t>;

	template <typename T, typename... Ts> concept same_as_any = (... or std::same_as<T, Ts>);


} // namespace april

