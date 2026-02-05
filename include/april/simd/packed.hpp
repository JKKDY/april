#pragma once
#include <april/simd/simd_traits.hpp>



#define AP_SIMD_BACKEND_XSIMD


//-----------------------
// SIMD BACKEND SWITCHING
//-----------------------
#if defined(AP_SIMD_BACKEND_XSIMD)
#include "april/simd/backend_xsimd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Packed = internal::xsimd::Packed<T, W>;
    template<typename T> using PackedMask = internal::xsimd::Mask<T>;
}

#elif  defined(AP_SIMD_BACKEND_STD_SIMD)

#include "april/simd/backend_std_simd.hpp"

namespace april::simd {
    template<typename T, size_t W = 0> using Packed = internal::std_simd::Packed<T, W>;
    template<typename T> using PackedMask = internal::std_simd::Mask<T>;
}
#endif


static_assert(april::simd::IsSimdType<april::simd::Packed<double>>);
static_assert(april::simd::IsSimdType<april::simd::Packed<float>>);
static_assert(april::simd::IsSimdMask<april::simd::PackedMask<double>>);
static_assert(april::simd::IsSimdMask<april::simd::PackedMask<float>>);
