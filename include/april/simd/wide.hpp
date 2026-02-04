#pragma once
#include <april/simd/concepts.hpp>




#if defined(AP_SIMD_BACKEND_XSIMD)
#include "april/simd/backend_xsimd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Wide = internal::xsimd::Wide<T, W>;
    template<typename T> using Mask = internal::xsimd::Mask<T>;
    AP_SIMD_IMPORT_WIDE_MATH(xsimd);
}

#elif  defined(AP_SIMD_BACKEND_STD_SIMD)

#include "april/simd/backend_std_simd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Wide = internal::std_simd::Wide<T, W>;
    template<typename T> using Mask = internal::std_simd::Mask<T>;
    AP_SIMD_IMPORT_WIDE_MATH(std_simd);
}

#endif


namespace april::simd {
    static_assert(IsSimdType<Wide<double>>);
    static_assert(IsSimdType<Wide<float>>);
}
