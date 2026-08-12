#pragma once
#include "april/base/backend_config.hpp"
#include "april/particle/properties.hpp"
#include "april/simd/packed_concept.hpp"



#if (defined(APRIL_SIMD_BACKEND_XSIMD) + \
defined(APRIL_SIMD_BACKEND_STD_SIMD) + \
defined(APRIL_SIMD_BACKEND_SCALAR) \
) > 1

#error "[APRIL] Multiple SIMD backends defined. Select exactly one."

#endif


//-----------------------
// SIMD BACKEND SWITCHING
//-----------------------
#if defined(APRIL_SIMD_BACKEND_XSIMD)

    #include "april/simd/backends/backend_xsimd.hpp"

    namespace april::simd {
        template<typename T, size_t W = 8>
        using Packed = internal::xsimd::Packed<T, W>;

        template<typename T, size_t W = 8>
        using PackedMask = internal::xsimd::Mask<T, W>;
    }


#elif defined(APRIL_SIMD_BACKEND_STD_SIMD)

    #include "april/simd/backends/backend_std_simd.hpp"

    namespace april::simd {
        template<typename T, size_t W = 2>
        using Packed = internal::std_simd::Packed<T, W>;

        template<typename T, size_t W = 2>
        using PackedMask = internal::std_simd::Mask<T, W>;
    }

#elif defined(APRIL_SIMD_BACKEND_SCALAR)

    #include "april/simd/backends/backend_scalar.hpp"

    namespace april::simd {
        template<typename T, size_t W = 8>
        using Packed = internal::scalar::Packed<T, W>;

        template<typename T, size_t W = 8>
        using PackedMask = internal::scalar::Mask<T, W>;
    }

#else
#error "[APRIL] No SIMD backend selected."
#endif

// static checks
namespace april::simd::internal {
    template<typename From, typename To>
    concept IsSimdMaskConvertibleTo =
        IsSimdMask<std::remove_cvref_t<From>> &&
        IsSimdMask<std::remove_cvref_t<To>> &&
        std::convertible_to<
            const std::remove_cvref_t<From>&,
            std::remove_cvref_t<To>
        >;
}

// need to import here so april::select is findable
#include "april/simd/math.hpp"

namespace april::simd::internal {
    template<typename Mask, typename Packed>
    concept IsSimdMaskCompatibleWith =
        IsSimdType<std::remove_cvref_t<Packed>> &&
        IsSimdMaskConvertibleTo<
               std::remove_cvref_t<Mask>,
            typename std::remove_cvref_t<Packed>::mask_type
        > &&
        requires(
               const std::remove_cvref_t<Mask>& mask,
            const std::remove_cvref_t<Packed>& a,
            const std::remove_cvref_t<Packed>& b
        ) {
            { april::select(mask, a, b) } -> std::same_as<std::remove_cvref_t<Packed>>;
        };

    static_assert(april::simd::IsSimdType<Packed<double>>);
    static_assert(april::simd::IsSimdType<Packed<float>>);
    static_assert(april::simd::IsSimdType<Packed<size_t>>);
    static_assert(april::simd::IsSimdType<Packed<int>>);
    static_assert(april::simd::IsSimdType<Packed<ParticleState>>);
    static_assert(april::simd::IsSimdType<Packed<ParticleID>>);
    static_assert(april::simd::IsSimdType<Packed<ParticleType>>);

    static_assert(april::simd::IsSimdMask<PackedMask<double>>);
    static_assert(april::simd::IsSimdMask<PackedMask<float>>);
    static_assert(april::simd::IsSimdMask<PackedMask<size_t>>);
    static_assert(april::simd::IsSimdMask<PackedMask<int>>);

    static_assert(IsSimdMaskConvertibleTo<PackedMask<float>, PackedMask<int>>);
    static_assert(IsSimdMaskConvertibleTo<PackedMask<int>, PackedMask<float>>);
    static_assert(IsSimdMaskConvertibleTo<PackedMask<double>, PackedMask<size_t>>);
    static_assert(IsSimdMaskConvertibleTo<PackedMask<size_t>, PackedMask<double>>);

    static_assert(IsSimdMaskCompatibleWith<PackedMask<int>, Packed<float>>);
    static_assert(IsSimdMaskCompatibleWith<PackedMask<float>, Packed<int>>);
    static_assert(IsSimdMaskCompatibleWith<PackedMask<size_t>, Packed<double>>);
    static_assert(IsSimdMaskCompatibleWith<PackedMask<double>, Packed<size_t>>);

    static_assert(std::same_as<Packed<float>::mask_type, PackedMask<float>>);
    static_assert(std::same_as<Packed<int>::mask_type, PackedMask<int>>);
    static_assert(std::same_as<Packed<double>::mask_type, PackedMask<double>>);
    static_assert(std::same_as<Packed<size_t>::mask_type, PackedMask<size_t>>);

    static_assert(IsSimdMaskConvertibleTo<PackedMask<float, 4>, PackedMask<int, 4>>);
    static_assert(IsSimdMaskCompatibleWith<PackedMask<int, 4>, Packed<float, 4>>);
    static_assert(!IsSimdMaskConvertibleTo<PackedMask<float, 4>, PackedMask<double, 2>>);
    static_assert(!IsSimdMaskCompatibleWith<PackedMask<double, 2>, Packed<float, 4>>);
}


