#pragma once

#include <concepts>
#include <type_traits>

#include "april/base/types.hpp"

// ---------------------------------------------------------
// TRIVIAL ATTRIBUTE MACRO
// ---------------------------------------------------------
// Stopgap macro to automatically generate the required C++20
// SIMD reflection boilerpate for single-field attributes.
//
// Usage:
// struct Charge {
//     APRIL_TRIVIAL_ATTRIBUTE(double /*type*/, charge /*name*/);
// };
// ---------------------------------------------------------
#define APRIL_TRIVIAL_ATTRIBUTE(ScalarT, FieldName) \
    ScalarT FieldName; \
    struct VectorLayout { \
        using ScalarType = ScalarT; \
        april::simd::Packed<ScalarT> FieldName; \
    }

namespace april {
    struct NoParticleAttributes {};
}

namespace april::particle {

    // ==== STOP GAP SOLUTION ====
    // this is a stopgap solution until reflection arrives
    // this allows Attributes with a single arithmetic type to be vectorized
    template <typename T>
     concept HasTrivialVectorLayout = requires {
        typename T::VectorLayout;
        typename T::VectorLayout::ScalarType;
     };

    template <typename T>
    concept IsTriviallyVectorizable = HasTrivialVectorLayout<T> && requires {
        // The underlying type must be arithmetic for SIMD math
        requires std::is_arithmetic_v<typename T::VectorLayout::ScalarType>;

        // The outer struct must be exactly the size of the scalar (e.g., 8 bytes for double)
        // This guarantees no hidden fields or padding
        requires sizeof(T) == sizeof(T::VectorLayout::ScalarType);

        // The layout struct must be exactly the size of one SIMD register
        // This guarantees they didn't add multiple packed registers to the layout
        requires sizeof(T::VectorLayout) == sizeof(simd::Packed<typename T::VectorLayout::ScalarType>);

        // Must be standard layout for safe reinterpret_cast
        requires std::is_standard_layout_v<T>;
        requires std::is_trivially_copyable_v<T>;
        requires std::is_standard_layout_v<typename T::VectorLayout>;
    };
    // ==== STOP GAP SOLUTION END ====

    template <typename T>
    concept IsVectorizable = IsTriviallyVectorizable<T>; // note: with C++26 reflection automatic best effort vectorization will be added


    template <typename T>
    concept IsParticleAttributes =
        std::default_initializable<T> &&
        std::is_trivially_copyable_v<T> &&
        std::is_trivially_destructible_v<T> &&
        std::is_standard_layout_v<T> &&
        (!std::is_polymorphic_v<T>);


    // template struct used to tell the environment what user data will be used
    template<typename Data = NoParticleAttributes>
    struct ParticleAttributes { using particle_attributes_t = Data; };
}

namespace april {
    template<typename Data = NoParticleAttributes>
     inline constexpr particle::ParticleAttributes<Data> particle_attributes {};
}

