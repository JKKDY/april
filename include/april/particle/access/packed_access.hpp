/**
 * @file packed_access.hpp
 * @brief SIMD (packed) particle access layer for vectorized kernels.
 *
 * This file mirrors scalar_access.hpp but operates on full SIMD registers instead of scalars.
 *
 * Key types and their roles:
 *
 * 1. PackedParticleRef     - "Vector pointer": locates a block of particles in AoSoA memory.
 *                            Operations on it usually cause immediate memory traffic.
 *
 * 2. PackedParticleBuffer  - Shadow copy in registers. This is where actual computation happens.
 *                            Allows multiple operations without touching memory repeatedly.
 *
 * 3. PackedBufferView      - Restricted view passed to the user kernel. Enforces that kernels
 *                            cannot write to read-only fields.
 *
 *
 * The separation between Ref / Buffer / View is deliberate: it maximizes register reuse
 * while keeping kernels simple and safe.
 */
#pragma once

#include "april/particle/access/internal/packed_buffer.hpp"
#include "april/particle/access/internal/packed_reference.hpp"
#include "april/particle/access/internal/packed_buffer_view.hpp"

namespace april::particle::internal {


    //---------
    // TRAITS
    //---------
    template <typename T>
    struct is_packed_buffer_impl : std::false_type {};

    template <typename T>
    struct is_packed_ref_impl : std::false_type {};

    template <typename T>
    struct is_buffer_view_impl : std::false_type {};

    // specialization for the unified types

    template <ParticleField RM, ParticleField WM, typename U>
    struct is_packed_ref_impl<PackedParticleRef<RM, WM, U>> : std::true_type {};

    template <typename Ref, typename Mask>
    struct is_packed_ref_impl<MaskedPackedParticleRef<Ref, Mask>> : std::true_type {};

    template <ParticleField RM, ParticleField WM, class Attributes>
    struct is_buffer_view_impl<PackedBufferView<RM, WM, Attributes>> : std::true_type {};
} // namespace april::particle::internal


namespace april::particle {
    //---------
    // CONCEPTS
    //---------
    /**
     * Concepts that identify the different kinds of particle accessors.
     *
     * These are the main interface contracts used throughout the library:
     * - When writing custom forces, boundaries, or controllers, you will usually
     *   see `IsAnyParticleAccessor` or `IsScalarParticleAccessor` in templates.
     */
    template <typename T>
    concept IsPackedParticleBuffer = internal::is_packed_buffer_impl<std::remove_cvref_t<T>>::value;

    template <typename T>
    concept IsPackedParticleRef = internal::is_packed_ref_impl<std::remove_cvref_t<T>>::value;

    template <typename T>
    concept IsPackedParticleView = internal::is_buffer_view_impl<std::remove_cvref_t<T>>::value;

    /**
     * Any packed (SIMD) accessor — buffer, ref, or view.
     * Most internal code uses this when it doesn't care about the exact form.
     */
    template <typename T>
    concept IsPackedParticleAccessor =
        IsPackedParticleBuffer<T> ||
        IsPackedParticleRef<T> ||
        IsPackedParticleView<T>;

    /**
     * Union of scalar and packed accessors.
     * This is the most commonly used concept when writing kernels.
     * A kernel can accept any accessor that provides p.position, p.velocity, etc.
     */
    template <typename T>
    concept IsAnyParticleAccessor = IsScalarParticleAccessor<T> || IsPackedParticleAccessor<T>;
} // namespace april::particle
