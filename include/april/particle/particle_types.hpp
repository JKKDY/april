#pragma once

#include <concepts>
#include <cstdint>
#include <type_traits>

#include "april/base/bitmask.hpp"

namespace april {

    /**
     * @brief Bitmask describing which particle fields a kernel or component may access.
     *
     * Field masks are used to expose only the requested particle data to kernels.
     * Accessing a field that was not declared in the corresponding mask is a
     * compile-time error.
     */
    enum class ParticleField : uint16_t {
        none         = 0u,
        position     = 1u << 0,
        velocity     = 1u << 1,
        force        = 1u << 2,
        old_position = 1u << 3,
        state        = 1u << 4,
        mass         = 1u << 5,
        type         = 1u << 6,
        id           = 1u << 7,
        attributes   = 1u << 8,

        all = position
            | velocity
            | force
            | old_position
            | state
            | mass
            | type
            | id
            | attributes
    };
    AP_ENABLE_BITMASK_OPERATORS(ParticleField)


    /**
     * @brief Bitmask describing particle activity during traversal and force evaluation.
     *
     * The base states describe how a particle participates in simulation.
     * Composite values such as EXERTING and MOVABLE are convenience filters.
     */
    enum class ParticleState : uint8_t {
        NONE       = 0u,
        ALIVE      = 1u << 0,  ///< Moves, exerts forces, and experiences forces.
        DEAD       = 1u << 1,  ///< Inactive; no movement and no interactions.
        PASSIVE    = 1u << 2,  ///< Moves and experiences forces, but exerts none.
        STATIONARY = 1u << 3,  ///< Exerts forces, but does not move or respond.
        INVALID    = 1u << 7,  ///< Sentinel for invalid storage slots or gaps.

        EXERTING = ALIVE | STATIONARY,
        MOVABLE = ALIVE | PASSIVE,

        ALL = ALIVE | DEAD | PASSIVE | STATIONARY
    };
    AP_ENABLE_BITMASK_OPERATORS(ParticleState)


    /// User-facing particle type used for interaction lookup.
    using ParticleType = uint16_t;

    /// User-facing particle ID used for stable particle references.
    using ParticleID = uint32_t;


    namespace particle::internal {

        template<class T>
        concept HasFields = requires {
            { std::remove_cvref_t<T>::fields } -> std::convertible_to<ParticleField>;
        };

        template<HasFields Self>
        inline constexpr ParticleField FieldOf = std::remove_cvref_t<Self>::fields;

        template<ParticleField M, ParticleField F>
        inline constexpr bool has_field_v = (M & F) != ParticleField::none;

    } // namespace particle::internal

} // namespace april