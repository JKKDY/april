#pragma once
#include "april/base/config.hpp"
#include "april/simd/simd_traits.hpp"



#if (defined(APRIL_SIMD_BACKEND_XSIMD) + \
defined(APRIL_SIMD_BACKEND_STD_SIMD)) > 1

#error "[APRIL] Multiple SIMD backends defined. Select exactly one."

#endif


//-----------------------
// SIMD BACKEND SWITCHING
//-----------------------
#if defined(APRIL_SIMD_BACKEND_XSIMD)

    #include "april/simd/backend_xsimd.hpp"

    namespace april::simd {
        template<typename T, size_t W = 0>
        using Packed = internal::xsimd::Packed<T, W>;

        template<typename T>
        using PackedMask = internal::xsimd::Mask<T>;
    }

#elif defined(APRIL_SIMD_BACKEND_STD_SIMD)

#include "april/simd/backend_std_simd.hpp"

namespace april::simd {
        template<typename T, size_t W = 0>
        using Packed = internal::std_simd::Packed<T, W>;

        template<typename T>
        using PackedMask = internal::std_simd::Mask<T>;
    }

#else
#error "[APRIL] No SIMD backend selected."
#endif


static_assert(april::simd::IsSimdType<april::simd::Packed<double>>);
static_assert(april::simd::IsSimdType<april::simd::Packed<float>>);
static_assert(april::simd::IsSimdType<april::simd::Packed<size_t>>);
static_assert(april::simd::IsSimdType<april::simd::Packed<int>>);

static_assert(april::simd::IsSimdMask<april::simd::PackedMask<double>>);
static_assert(april::simd::IsSimdMask<april::simd::PackedMask<float>>);
static_assert(april::simd::IsSimdMask<april::simd::PackedMask<size_t>>);
static_assert(april::simd::IsSimdMask<april::simd::PackedMask<int>>);


