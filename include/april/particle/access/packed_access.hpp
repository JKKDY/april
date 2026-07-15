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
#include "april/base/types.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/access/source.hpp"
#include "april/particle/properties.hpp"
#include "scalar_access.hpp"
#include "april/particle/attributes.hpp"

namespace april::particle::internal {
    // forward declaration
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedBufferView;
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedParticleRef;
    template <typename Ref, typename Mask> struct MaskedPackedParticleRef;

    //--------------------
    // PACKED PARTICLE REF
    //--------------------
    /**
     * SIMD equivalent of ScalarParticleRef.
     *
     * Holds packed (SIMD-width) pointers/references to a contiguous block of particles
     * in AoSoA layout. It does *not* load data into registers yet. That happens in load_buffer().
     *
     * This struct is intentionally lightweight so the compiler can inline it aggressively.
     */
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
    struct PackedParticleRef {
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask & ~ParticleField::id;

    private:
        /**
          * Initializes a packed field pointer from the ParticleSource.
          * Special handling for enums (state, type) to convert them to their underlying integer type
          * while preserving const-correctness for SIMD compatibility.
          */
        template <ParticleField F, typename Source>
        constexpr auto init_packed(const Source& src) {
            if constexpr (particle::internal::has_field_v<ReadAccess | WriteAccess, F>) {
                auto ptr = src.template get<F>();

                // get value type from pointer
                using PtrType = decltype(ptr);
                using ValueType = std::remove_pointer_t<PtrType>;

                if constexpr (std::is_enum_v<ValueType>) {
                    // Determine the underlying integer and preserve constness
                    using IntType = std::underlying_type_t<ValueType>;
                    using TargetPtr = std::conditional_t<std::is_const_v<ValueType>, const IntType*, IntType*>;

                    return reinterpret_cast<TargetPtr>(ptr);
                } else {
                    // Return normal pointers as-is
                    return ptr;
                }
            }
            else {
                return AccessForbidden<F>();
            }
        }

        // Resolves the member type to Mutable, Const, or Poison based on masks
        template <typename MutT, typename ConstT, ParticleField F>
        using field_t = field_access_t<MutT, ConstT, F, ReadAccess, WriteAccess>;

        // Maps scalars (double, int) to their corresponding SIMD register types
        template <typename T>
        using target_reg_t = std::conditional_t<
            std::is_floating_point_v<T>,
            packed::value_type,
            std::conditional_t<std::is_signed_v<T>, packedi::value_type, packedu::value_type>
        >;

        // Helper for single-component packed fields (mass, type, etc.)
        template <typename MemT, ParticleField F>
        using packed_field_t = field_t<
            simd::PackedRef<MemT, simd::Packed<target_reg_t<MemT>>>,
            const simd::PackedRef<const MemT, simd::Packed<target_reg_t<MemT>>>,
            F
        >;

        // Helper for 3D vector packed fields (pos, vel, etc.)
        template <ParticleField F>
        using vec3_field_t = field_t<math::Vec3Proxy<pvec3::type>, const math::Vec3Proxy<const pvec3::type>, F>;

        // Declare a raw pointer
        template<typename T, ParticleField F>
        using Ptr = field_access_t<T* APRIL_RESTRICT, const T* APRIL_RESTRICT, F, ReadAccess, WriteAccess>;

    public:
        explicit PackedParticleRef(const auto& source) noexcept
           : force(init_packed<ParticleField::force>(source))
           , position(init_packed<ParticleField::position>(source))
           , velocity(init_packed<ParticleField::velocity>(source))
           , old_position(init_packed<ParticleField::old_position>(source))
           , mass(init_packed<ParticleField::mass>(source))
           , state(init_packed<ParticleField::state>(source))
           , type(init_packed<ParticleField::type>(source))
           , id(init_packed<ParticleField::id>(source))
           , attributes(init_packed<ParticleField::attributes>(source))
            {}

        /**
         * Narrowing constructor — used to create read-only views from mutable references.
         * Expansion of write permissions is forbidden at compile time.
         */
        template <ParticleField OtherWriteMask>
            requires ((WriteAccess & OtherWriteMask) == WriteAccess)
        explicit PackedParticleRef(const PackedParticleRef<ReadAccess, OtherWriteMask, Attributes>& r) noexcept
            : force(r.force)
              , position(r.position)
              , velocity(r.velocity)
              , old_position(r.old_position)
              , mass(r.mass)
              , state(r.state)
              , type(r.type)
              , id(r.id)
              , attributes(r.attributes)
        {}

        /**
          * Returns a read-only view of this SIMD block.
          */
        auto to_view() const noexcept {
            return PackedParticleRef<ReadAccess | WriteAccess, ParticleField::none, Attributes>(*this);
        }

        /**
         * Loads the referenced memory block into SIMD registers (PackedParticleBuffer).
         * This is the main transition point from memory to computation.
         */
        PackedParticleBuffer<ReadAccess, WriteAccess, Attributes> load_buffer() const noexcept {
            return PackedParticleBuffer<ReadAccess, WriteAccess, Attributes>(*this);
        }

        /**
          * Future stub for applying persistent SIMD masks.
          */
        template<typename Mask>
        auto mask_with(const Mask & mask) {
            return MaskedPackedParticleRef<std::remove_cvref_t<decltype(*this)>, Mask>(*this, mask);
        }

        // Data members with strict const-correctness
        // APRIL_NO_UNIQUE_ADDRESS ensures forbidden fields do not increase object size.
        APRIL_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::force> force;
        APRIL_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::position> position;
        APRIL_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::velocity> velocity;
        APRIL_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::old_position> old_position;

        APRIL_NO_UNIQUE_ADDRESS packed_field_t<double, ParticleField::mass> mass;
        APRIL_NO_UNIQUE_ADDRESS packed_field_t<std::underlying_type_t<ParticleState>, ParticleField::state> state;
        APRIL_NO_UNIQUE_ADDRESS packed_field_t<ParticleType, ParticleField::type> type;
        APRIL_NO_UNIQUE_ADDRESS packed_field_t<ParticleID, ParticleField::id> id;

        APRIL_NO_UNIQUE_ADDRESS Ptr<Attributes, ParticleField::attributes> attributes;
    };

    /**
     * @brief Decorator for propagating SIMD masks through layered calls.
     * * This allows filters (e.g., state checks) to be applied once and
     * carried down through the kernel call chain, avoiding redundant
     * mask re-computation in deeply nested user code.
     * Currently unused as it is a future stub.
     */
    template <typename Ref, typename Mask>
    struct MaskedPackedParticleRef : Ref {
        Mask mask;

        explicit MaskedPackedParticleRef(const Ref& r, const Mask& m) noexcept
            : Ref(r), mask(m) {}

        auto mask_with(const Mask & m) const noexcept {
            return MaskedPackedParticleRef{static_cast<const Ref&>(*this), mask & m};
        }
    };



    //-----------------------
    // PACKED PARTICLE BUFFER
    //-----------------------
    /**
     * Register-backed "shadow" object for a block of particles.
     *
     * DESIGN INTENT:
     * We load particle data once into SIMD registers, perform all interactions
     * entirely in registers (to minimize memory traffic), and write back exactly once.
     *
     * MASKING STRATEGY:
     * - RWMask (Read+Write): Fields that are overwritten (e.g. position, velocity).
     * - WOMask (Write-Only): Fields that accumulate deltas (e.g. force).
     *   These are zero-initialized so they act as clean accumulators.
     * - ROMask (Read-Only): Constant fields (e.g. mass, type).
     */
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
    struct PackedParticleBuffer {
    private:
        // if in ReadMask or WriteMask return type T else return Poison
        template <typename T, ParticleField F>
          using buffer_field_t = std::conditional_t<
              has_field_v<ReadMask | WriteMask, F>,
              T,
              AccessForbidden<F>
          >;

        template <ParticleField F>
        using pvec3_t = buffer_field_t<pvec3, F>;

        template <ParticleField F>
        using packed_float_t = buffer_field_t<packed, F>;

        template <ParticleField F>
        using packed_int_t = buffer_field_t<packedu, F>;

        // ==== STOP GAP SOLUTION ==== (will be replaced in C++26 with reflection)
        // vectorization of simple scalar type
        static_assert(!has_field_v<ReadMask | WriteMask, ParticleField::attributes> || IsTriviallyVectorizable<Attributes>,
            "Vectorization of non trivial attributes not possible yet");
        // traits case for scalar extraction
        template <typename T>
        struct extract_attr_scalar { using type = double; /*Fallback*/ };

        template <typename T> requires IsTriviallyVectorizable<T>
        struct extract_attr_scalar<T> { using type = T::VectorLayout::ScalarType; };

        using attr_scalar_t = extract_attr_scalar<Attributes>::type;

        template <ParticleField F>
        using packed_attr_t = buffer_field_t<simd::Packed<attr_scalar_t>, F>;
        // ====== STOP GAP SOLUTION END =====

    public:
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask & ~ParticleField::id; // soft restriction so user does not need to zero out the id specifically when using ParticleField::all

        static constexpr ParticleField RWMask = ReadMask & WriteMask;  // read & write
        static constexpr ParticleField WOMask = WriteMask & ~ReadMask; // write only
        static constexpr ParticleField ROMask = ReadMask & ~WriteMask; // read only

        APRIL_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::position> position;
        APRIL_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::old_position> old_position;
        APRIL_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::velocity> velocity;
        APRIL_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::force> force;

        APRIL_NO_UNIQUE_ADDRESS packed_float_t<ParticleField::mass> mass;
        APRIL_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::state> state;
        APRIL_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::type>  type;
        APRIL_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::id>    id;

        APRIL_NO_UNIQUE_ADDRESS packed_attr_t<ParticleField::attributes> attributes;


        PackedParticleBuffer() = default;

        /**
         * Load from memory (PackedParticleRef).
         *
         * ReadMask fields are loaded from memory.
         * WOMask fields are zero-initialized to serve as clean accumulators.
         * This is critical for symmetric interactions using the rotation sweep,
         * where each force contribution is added reciprocally (see container/batching/chunked_batch.hpp).
         */
        template <typename attr>
        explicit PackedParticleBuffer(const PackedParticleRef<ReadMask, WriteMask, attr>& source) {
            // Load Read-enabled fields
            if constexpr (has_field_v<ReadMask, ParticleField::position>) position = source.position;
            if constexpr (has_field_v<ReadMask, ParticleField::old_position>) old_position = source.old_position;
            if constexpr (has_field_v<ReadMask, ParticleField::velocity>) velocity = source.velocity;
            if constexpr (has_field_v<ReadMask, ParticleField::force>) force = source.force;
            if constexpr (has_field_v<ReadMask, ParticleField::mass>) mass = source.mass;
            if constexpr (has_field_v<ReadMask, ParticleField::state>) state = source.state;
            if constexpr (has_field_v<ReadMask, ParticleField::type>) type = source.type;
            if constexpr (has_field_v<ReadMask, ParticleField::id>) id = source.id;

            if constexpr (has_field_v<ReadMask, ParticleField::attributes>) {
                // Cast the AoS struct pointer to an arithmetic pointer
                // IsTriviallyVectorizable ensures that this reinterpret casting will work fine
                auto ptr = reinterpret_cast<const attr_scalar_t*>(source.attributes);
                attributes = decltype(attributes)::load_unaligned(ptr);
            }

            // Zero-initialize Write-Only accumulators (WOMask).
            // This transforms the register into a "delta-buffer" for numerical types
            if constexpr (has_field_v<WOMask, ParticleField::position>) position = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) old_position = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) velocity = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::force>) force = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::mass>) mass = 0.0;

            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                attributes = decltype(attributes)(0);
            }
        }

        /**
         * Broadcast a single scalar particle into all lanes.
         * Used for 1 × N interactions (e.g. tail particle vs full block).
         */
        template <typename ScalarAccessor>
            requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        explicit PackedParticleBuffer(const ScalarAccessor& scalar) {
            auto broad_cast_vec = [&]<ParticleField field>(auto && packed_vec, auto && scalar_vec) APRIL_FORCE_INLINE {
                if constexpr (has_field_v<ReadMask, field>) {
                    packed_vec.x = scalar_vec.x;
                    packed_vec.y = scalar_vec.y;
                    packed_vec.z = scalar_vec.z;
                } else if constexpr (has_field_v<WOMask, field>) {
                    packed_vec = pvec3(0.0);
                }
            };

            broad_cast_vec.template operator()<ParticleField::position>(position, scalar.position);
            broad_cast_vec.template operator()<ParticleField::old_position>(old_position, scalar.old_position);
            broad_cast_vec.template operator()<ParticleField::velocity>(velocity, scalar.velocity);
            broad_cast_vec.template operator()<ParticleField::force>(force, scalar.force);

            if constexpr (has_field_v<ReadMask, ParticleField::mass>) {
                mass = scalar.mass;
            } else if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass = 0.0;
            }

            if constexpr (has_field_v<ReadMask, ParticleField::state>) state = scalar.state;
            if constexpr (has_field_v<ReadMask, ParticleField::type>) type = scalar.type;
            if constexpr (has_field_v<ReadMask, ParticleField::id>) id = scalar.id;

            if constexpr (has_field_v<ReadMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<const attr_scalar_t*>(&scalar.attributes);
                attributes = decltype(attributes)(*ptr);
            }
        }

        // export as view
        APRIL_FORCE_INLINE PackedBufferView<ReadMask, WriteMask, Attributes> to_view() {
            return PackedBufferView(*this);
        }


        /**
          * Register Lane Rotations.
          */
        template <unsigned K = 1>
        APRIL_FORCE_INLINE void rotate_left() {
            auto rotate_vec = [&]<ParticleField field>(auto && vec) APRIL_FORCE_INLINE {
                if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, field>) {
                    vec.x = vec.x.template rotate_left<K>();
                    vec.y = vec.y.template rotate_left<K>();
                    vec.z = vec.z.template rotate_left<K>();
                }
            };

            auto rotate_scalar = [&]<ParticleField field>(auto && scalar) APRIL_FORCE_INLINE {
                if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, field>) {
                    scalar = scalar.template rotate_left<K>();
                }
            };

            rotate_vec.template operator()<ParticleField::position>(position);
            rotate_vec.template operator()<ParticleField::old_position>(old_position);
            rotate_vec.template operator()<ParticleField::velocity>(velocity);
            rotate_vec.template operator()<ParticleField::force>(force);

            rotate_scalar.template operator()<ParticleField::mass>(mass);
            rotate_scalar.template operator()<ParticleField::state>(state);
            rotate_scalar.template operator()<ParticleField::type>(type);
            rotate_scalar.template operator()<ParticleField::id>(id);
            rotate_scalar.template operator()<ParticleField::attributes>(attributes);
        }
        template <unsigned K = 1>
        APRIL_FORCE_INLINE void rotate_right() {
            auto rotate_vec = [&]<ParticleField field>(auto && vec) APRIL_FORCE_INLINE {
                if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, field>) {
                    vec.x = vec.x.template rotate_right<K>();
                    vec.y = vec.y.template rotate_right<K>();
                    vec.z = vec.z.template rotate_right<K>();
                }
            };

            auto rotate_scalar = [&]<ParticleField field>(auto && scalar) APRIL_FORCE_INLINE {
                if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, field>) {
                    scalar = scalar.template rotate_right<K>();
                }
            };

            rotate_vec.template operator()<ParticleField::position>(position);
            rotate_vec.template operator()<ParticleField::old_position>(old_position);
            rotate_vec.template operator()<ParticleField::velocity>(velocity);
            rotate_vec.template operator()<ParticleField::force>(force);

            rotate_scalar.template operator()<ParticleField::mass>(mass);
            rotate_scalar.template operator()<ParticleField::state>(state);
            rotate_scalar.template operator()<ParticleField::type>(type);
            rotate_scalar.template operator()<ParticleField::id>(id);
            rotate_scalar.template operator()<ParticleField::attributes>(attributes);
        }

        /**
         * Accumulate reciprocal deltas from another buffer.
         * Only affects WOMask fields. Base state (RW/RO) is never overwritten.
         */
        APRIL_FORCE_INLINE void accumulate(const PackedParticleBuffer& other) {
            if constexpr (has_field_v<WOMask, ParticleField::position>) position += other.position;
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) old_position += other.old_position;
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) velocity += other.velocity;
            if constexpr (has_field_v<WOMask, ParticleField::force>) force += other.force;
            if constexpr (has_field_v<WOMask, ParticleField::mass>) mass += other.mass;
            if constexpr (has_field_v<WOMask, ParticleField::attributes>) attributes += other.attributes;
        }

        /**
        * Masked Delta Accumulation.
        */
        template <typename MaskT>
        APRIL_FORCE_INLINE void accumulate(const PackedParticleBuffer& other, const MaskT& mask) {
            const packed null = 0.0;

            auto accumulate_vec = [&]<ParticleField field>(auto && this_field, auto && other_field) APRIL_FORCE_INLINE {
                if constexpr (has_field_v<WOMask, field>) {
                    this_field.x += select(mask, other_field.x, packed(0));
                    this_field.y += select(mask, other_field.y, packed(0));
                    this_field.z += select(mask, other_field.z, packed(0));
                }
            };

            accumulate_vec.template operator()<ParticleField::position>(position, other.position);
            accumulate_vec.template operator()<ParticleField::old_position>(old_position, other.old_position);
            accumulate_vec.template operator()<ParticleField::velocity>(velocity, other.velocity);
            accumulate_vec.template operator()<ParticleField::force>(force, other.force);

            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass += select(mask, other.mass, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                attributes += select(mask, other.attributes, decltype(attributes)(0));
            }
        }

        /**
         * Flush register values back to memory.
         *
         * - WOMask fields: additive update (dest += src)
         * - RWMask fields: replacement update (dest = src)
         *
         * This distinction allows both in-place modification and force accumulation
         * to work correctly within the same framework.
         */
        template <typename Attr>
        APRIL_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr>& packed_ref) const {
            // Write-Only fields use additive accumulation (preserves base state)

            auto update_field = [&]<ParticleField Field>(auto&& dest, auto&& src) {
                if constexpr (has_field_v<WOMask, Field>) {
                    dest += src;
                } else if constexpr (has_field_v<RWMask, Field>) {
                    dest = src;
                }
            };

            update_field.template operator()<ParticleField::position>(packed_ref.position, position);
            update_field.template operator()<ParticleField::old_position>(packed_ref.old_position, old_position);
            update_field.template operator()<ParticleField::velocity>(packed_ref.velocity, velocity);
            update_field.template operator()<ParticleField::force>(packed_ref.force, force);
            update_field.template operator()<ParticleField::mass>(packed_ref.mass, mass);

            if constexpr (has_field_v<RWMask, ParticleField::state>) packed_ref.state = state;
            if constexpr (has_field_v<RWMask, ParticleField::type>) packed_ref.type = type;

            // ATTRIBUTES
            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(packed_ref.attributes);
                auto current =  decltype(attributes)::load_unaligned(ptr);
                (current + attributes).store_unaligned(ptr);
            } else if constexpr (has_field_v<RWMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(packed_ref.attributes);
                attributes.store_unaligned(ptr);
            }

            // id is not assignable
        }

        /**
         * Masked Memory Flush.
         */
        template <typename Attr, typename MaskT>
        APRIL_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr>& packed_ref, const MaskT & mask) const {
            update_vec_masked<ParticleField::position>(packed_ref.position, position, mask);
            update_vec_masked<ParticleField::old_position>(packed_ref.old_position, old_position, mask);
            update_vec_masked<ParticleField::velocity>(packed_ref.velocity, velocity, mask);
            update_vec_masked<ParticleField::force>(packed_ref.force, force, mask);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                packed_ref.mass += select(mask, mass, 0.0);
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                packed_ref.mass = select(mask, mass, packed_ref.mass);

            // ATTRIBUTES (Masked read-modify-write)
            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(packed_ref.attributes);
                auto current = decltype(attributes)::load_unaligned(ptr);
                auto updated = current + select(mask, attributes, 0.0);
                updated.store_unaligned(ptr);
            } else if constexpr (has_field_v<RWMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(packed_ref.attributes);
                auto current = decltype(attributes)::load_unaligned(ptr);
                auto updated = select(mask, attributes, current);
                updated.store_unaligned(ptr);
            }
        }

        // future stub
        template <typename Ref, typename Mask>
        APRIL_FORCE_INLINE void update_into(MaskedPackedParticleRef<Ref, Mask>& masked_ref) const {
            update_into(masked_ref, masked_ref.mask);
        }

        /**
         * Horizontal Reduction.
         * * Collapses all lanes in the SIMD register into a single scalar value/vector.
         * Used to finalize 1 x N interactions where forces from an entire chunk
         * are reduced into a single scalar "tail" particle.
         */
        template <typename ScalarAccessor>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        APRIL_FORCE_INLINE void reduce_into(ScalarAccessor& p) const {
            reduce_vec_unmasked<ParticleField::position>(p.position, position);
            reduce_vec_unmasked<ParticleField::old_position>(p.old_position, old_position);
            reduce_vec_unmasked<ParticleField::velocity>(p.velocity, velocity);
            reduce_vec_unmasked<ParticleField::force>(p.force, force);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                p.mass += mass.reduce_add();
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot reduce RW mass.");

            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(&p.attributes);
                *ptr += attributes.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::attributes>) {
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot reduce RW attributes.");
            }
        }

        /**
        * Masked Horizontal Reduction.
        */
        template <typename ScalarAccessor, typename MaskT>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        APRIL_FORCE_INLINE void reduce_into(ScalarAccessor& p, const MaskT& mask) const {
            reduce_vec_masked<ParticleField::position>(p.position, position, mask);
            reduce_vec_masked<ParticleField::old_position>(p.old_position, old_position, mask);
            reduce_vec_masked<ParticleField::velocity>(p.velocity, velocity, mask);
            reduce_vec_masked<ParticleField::force>(p.force, force, mask);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                p.mass += select(mask, mass, 0.0).reduce_add();
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot masked reduce RW mass.");

            // ATTRIBUTES
            if constexpr (has_field_v<WOMask, ParticleField::attributes>) {
                auto ptr = reinterpret_cast<attr_scalar_t*>(&p.attributes);
                *ptr += select(mask, attributes, decltype(attributes)(0)).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::attributes>) {
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot masked reduce RW attributes.");
            }
        }

    private:
        // Unified vector write-back (masked) for WO and WR fields
        template <ParticleField F, typename DestT, typename SrcT, typename MaskT>
        APRIL_FORCE_INLINE static void update_vec_masked(DestT& dest, const SrcT& src, const MaskT& mask) {
            if constexpr (has_field_v<WOMask, F>) {
                const packed null = 0.0;
                dest.x += select(mask, src.x, null);
                dest.y += select(mask, src.y, null);
                dest.z += select(mask, src.z, null);
            } else if constexpr (has_field_v<RWMask, F>) {
                dest.x = select(mask, src.x, dest.x);
                dest.y = select(mask, src.y, dest.y);
                dest.z = select(mask, src.z, dest.z);
            }
        }

        // Unified vector reduce (unmasked) for WO and WR fields
        template <ParticleField F, typename ScalarT, typename SimdT>
        APRIL_FORCE_INLINE static void reduce_vec_unmasked(ScalarT& dest, const SimdT& src) {
            if constexpr (has_field_v<WOMask, F>) {
                dest.x += src.x.reduce_add();
                dest.y += src.y.reduce_add();
                dest.z += src.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, F>) {
                static_assert(sizeof(ScalarT) == 0, "FATAL: Cannot reduce a Read-Write vector field from a SIMD register to a scalar.");
            }
        }

        // Unified vector reduce (masked) for WO and WR fields
        template <ParticleField F, typename ScalarT, typename SimdT, typename MaskT>
        APRIL_FORCE_INLINE static void reduce_vec_masked(ScalarT& dest, const SimdT& src, const MaskT& mask) {
            if constexpr (has_field_v<WOMask, F>) {
                const packed null = 0.0;
                dest.x += select(mask, src.x, null).reduce_add();
                dest.y += select(mask, src.y, null).reduce_add();
                dest.z += select(mask, src.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, F>) {
                static_assert(sizeof(ScalarT) == 0, "FATAL: Cannot perform masked reduction on a Read-Write vector field.");
            }
        }
    };


    //------------
    // BUFFER VIEW
    //------------
    /**
     * Restricted view of PackedParticleBuffer passed to user kernels.
     *
     * This is the final safety layer in the SIMD particle access system.
     * It transforms the raw register data from PackedParticleBuffer into:
     *   - Mutable references (T&)   for fields in WriteMask
     *   - Const references (const T&) for fields that are read-only
     *   - AccessForbidden poison   for any field not explicitly requested
     *
     * DESIGN GOAL:
     * Give kernels a natural syntax (`p.position`, `p.force += ...`) while
     * letting the compiler enforce the declared Read/Write contract at compile time.
     * Any illegal write will simply fail to compile with a clear error.
     *
     * Because PackedBufferView is essentially a thin bundle of references,
     * it has zero runtime overhead and is typically completely optimized away.
     */
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
    struct PackedBufferView {
    private:
        /**
         * Chooses the correct reference type for each field based on the WriteMask:
         *   - Mutable reference if writable
         *   - Const reference if read-only
         *   - Poison type if forbidden
         */
        template <ParticleField F, typename T>
          using view_ref_t = std::conditional_t<
              std::is_same_v<T, AccessForbidden<F>>, // if it's poison, keep it as poison (by value)
              T,
              std::conditional_t<has_field_v<WriteMask, F>, T&, const T&>
              // if it's valid and in WriteMask, make it a mutable ref else a const ref
          >;

        using Buffer = PackedParticleBuffer<ReadMask, WriteMask, Attributes>;


        // ==== STOP GAP SOLUTION ==== (will be replaced in C++26 with reflection)
        // Attribute Vectorization Logic of trivial scalar types
        template <typename T>
        struct extract_attr_vector { using type = void; /*Fallback*/ };

        template <typename T> requires IsTriviallyVectorizable<T>
        struct extract_attr_vector<T> { using type = T::VectorLayout; };

        using attr_vector_t = extract_attr_vector<Attributes>::type;

        // resolve final exposed type based on the Buffer's actual permissions
        using exposed_attr_t = std::conditional_t<
            std::is_same_v<decltype(Buffer::attributes), AccessForbidden<ParticleField::attributes>>,
            AccessForbidden<ParticleField::attributes>,
            attr_vector_t
        >;

        static decltype(auto) bind_attributes(Buffer& buf) {
            if constexpr (std::is_same_v<exposed_attr_t, AccessForbidden<ParticleField::attributes>>) {
                return AccessForbidden<ParticleField::attributes>{};
            } else {
                // Pointer-Interconvertibility: Cast the raw Packed<T> to VectorLayout
                // Safe reinterpretation thanks to standard layout guarantees
                return reinterpret_cast<exposed_attr_t&>(buf.attributes);
            }
        }
        // ==== STOP GAP SOLUTION END ====

    public:
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask & ~ParticleField::id;

        // Mapping of Buffer registers to View references.
		// APRIL_NO_UNIQUE_ADDRESS ensures forbidden fields do not increase object size.
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::position, decltype(Buffer::position)> position;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::old_position, decltype(Buffer::old_position)> old_position;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::velocity, decltype(Buffer::velocity)> velocity;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::force, decltype(Buffer::force)> force;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::mass, decltype(Buffer::mass)> mass;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::state, decltype(Buffer::state)> state;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::type, decltype(Buffer::type)> type;
        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::id, decltype(Buffer::id)> id;

        APRIL_NO_UNIQUE_ADDRESS view_ref_t<ParticleField::attributes, exposed_attr_t> attributes;

        /**
          * Binds the buffer's registers to the view's references.
          * Extremely lightweight — usually completely elided by the optimizer.
          */
        APRIL_FORCE_INLINE explicit PackedBufferView(Buffer& buf)
            : position(buf.position),
              old_position(buf.old_position),
              velocity(buf.velocity),
              force(buf.force),
              mass(buf.mass),
              state(buf.state),
              type(buf.type),
              id(buf.id),
              attributes(bind_attributes(buf))
            {}
    };



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
    template <ParticleField RM, ParticleField WM, class Attributes>
    struct is_packed_buffer_impl<PackedParticleBuffer<RM, WM, Attributes>> : std::true_type {};

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
