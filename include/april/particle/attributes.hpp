#pragma once

#include <concepts>
#include <type_traits>

#include "april/base/types.hpp"

namespace april {
    struct NoParticleAttributes {
        // empty layout
        struct VectorLayout {
            struct Registers {};
            static Registers load(const NoParticleAttributes*) { return {}; }
            static void store(NoParticleAttributes*, const Registers&, const auto&) {}
        };
    };
}

namespace april::particle {
    namespace internal {
        // resolve enums to their underlying arithmetic type, or leave as-is
        template <typename T>
        using underlying_arithmetic_t = std::conditional_t<std::is_enum_v<T>, std::underlying_type_t<T>, T>;
    }

    // 1. Trivial Vectorization Concept
    // Matches standard layout structs containing exactly one arithmetic type (or enum)
    template <typename T>
    concept IsTriviallyVectorizable = requires {
        typename T::VectorLayout;
    } &&
    std::is_arithmetic_v<internal::underlying_arithmetic_t<typename T::VectorLayout>> &&
    std::is_standard_layout_v<T> &&
    std::is_trivially_copyable_v<T> &&
    sizeof(T) == sizeof(T::VectorLayout);


    template <typename Layout, typename AttrT>
    concept IsVectorLayout = requires(
        const AttrT* cptr, AttrT* ptr, typename Layout::Registers regs, const typename Layout::Registers cregs, packed_mask m
    ) {
        // Must define the bundle of registers (e.g., a struct of Packed types)
        typename Layout::Registers;

        // Must provide a way to load from the AoS array into the Register bundle
        { Layout::load(cptr) } -> std::same_as<typename Layout::Registers>;

        // Must provide a way to store the bundle back to AoS memory
        // Note: We use the mask to ensure only active particles are updated
        { Layout::store(ptr, cregs, m) } -> std::same_as<void>;
    };

    // Check if the VectorLayout struct even exists
    template <typename T>
    concept HasVectorLayout = requires {
        typename T::VectorLayout;
    } && IsVectorLayout<typename T::VectorLayout, T>;

    template <typename T>
    concept IsVectorizable = IsTriviallyVectorizable<T> || HasVectorLayout<T>;


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

