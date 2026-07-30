#pragma once

#include "april/base/types.hpp"
#include "april/simd/packed.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/access/source.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/access/scalar_access.hpp"
#include "april/particle/attributes.hpp"


namespace april::particle::internal {
    // forward declaration
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedBufferView;
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedParticleRef;
    template <typename Ref, typename Mask> struct MaskedPackedParticleRef;


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

            rotate(rotate_vec, rotate_scalar);
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

            rotate(rotate_vec, rotate_scalar);
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
         * This distinction allows both in-place modification and accumulation
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
                static_assert(sizeof(ScalarT) == 0, "Cannot reduce a Read-Write vector field from a SIMD register to a scalar.");
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
                static_assert(sizeof(ScalarT) == 0, "Cannot perform masked reduction on a Read-Write vector field.");
            }
        }


        template<typename RotateVec, typename RotateScalar>
        APRIL_FORCE_INLINE void rotate(RotateVec && rotate_vec, RotateScalar && rotate_scalar) {
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


    };
}