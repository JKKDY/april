#pragma once

#include "april/base/types.hpp"
#include "april/simd/packed.hpp"
#include "april/particle/access/source.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"



namespace april::particle::internal {
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedParticleRef;


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

        packed_mask & write_mask;
        bool & is_masked;

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
        APRIL_FORCE_INLINE explicit PackedBufferView(Buffer& buf, packed_mask & write_mask, bool & is_masked):
            write_mask(write_mask),
            is_masked(is_masked),
            position(buf.position),
            old_position(buf.old_position),
            velocity(buf.velocity),
            force(buf.force),
            mass(buf.mass),
            state(buf.state),
            type(buf.type),
            id(buf.id),
            attributes(bind_attributes(buf))
        {}

        void mask_with(const packed_mask new_mask) const {
            if (is_masked) {
                write_mask = write_mask & new_mask;
            } else {
                write_mask = new_mask;
                is_masked = true;
            }
        }
    };

}

