#pragma once

#include "april/base/types.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
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

}
