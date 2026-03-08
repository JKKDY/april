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
        // 1. The underlying type must be arithmetic
        requires std::is_arithmetic_v<typename T::VectorLayout::ScalarType>;

        // 2. The outer struct (Scalar/AoS) must be exactly the size of one scalar
        // This is the 8-byte check for a double.
        requires sizeof(T) == sizeof(typename T::VectorLayout::ScalarType);

        // 3. The inner VectorLayout struct must be exactly the size of one SIMD register
        // This is the 32-byte check for AVX.
        requires sizeof(typename T::VectorLayout) == sizeof(typename simd::Packed<typename T::VectorLayout::ScalarType>);

        // 4. Standard layout requirements for pointer interconvertibility
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

