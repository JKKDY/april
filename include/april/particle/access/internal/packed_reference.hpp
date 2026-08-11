#pragma once

#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/access/source.hpp"
#include "april/particle/properties.hpp"
#include "april/particle/attributes.hpp"
#include "april/particle/access/internal/fwd.hpp"

namespace april::particle::internal {


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
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes, typename Source>
    struct PackedParticleRef {
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask & ~ParticleField::id;

    private:
        template<ParticleField F>
        using raw_source_field_t = std::remove_cvref_t<decltype(std::declval<const Source&>().template get<F>())>;

        template<ParticleField F>
        static constexpr auto init_field(const Source & source) {
            if constexpr (!has_field_v<ReadAccess | WriteAccess, F>) {
                return AccessForbidden<F>();
            } else {
                auto source_field = source.template get<F>();
                using source_field_t = std::remove_cvref_t<decltype(source_field)>;

                auto to_location = []<typename FieldT>(const FieldT & field) {
                    if constexpr (std::is_pointer_v<FieldT>) {
                        return  simd::ContiguousLocation<std::remove_pointer_t<FieldT>>(field);
                    } else if constexpr (simd::IsLocation<FieldT>){
                        return field;
                    } else {
                        static_assert(false, "field type is not a raw pointer, nor a valid Location");
                    }
                };

                if constexpr (math::IsVec3Location<source_field_t>){
                    return math::Vec3Proxy(
                        simd::PackedRef(to_location(source_field.x)),
                        simd::PackedRef(to_location(source_field.y)),
                        simd::PackedRef(to_location(source_field.z))
                    );
                } else {
                    return  simd::PackedRef(to_location(source_field));
                }
            }
        }

        template<ParticleField F>
        using field_t = decltype(init_field<F>(std::declval<const Source&>()));

    public:
        explicit PackedParticleRef(const Source& source) noexcept
          : force(init_field<ParticleField::force>(source))
          , position(init_field<ParticleField::position>(source))
          , velocity(init_field<ParticleField::velocity>(source))
          , old_position(init_field<ParticleField::old_position>(source))
          , mass(init_field<ParticleField::mass>(source))
          , state(init_field<ParticleField::state>(source))
          , type(init_field<ParticleField::type>(source))
          , id(init_field<ParticleField::id>(source))
        {}

        /**
         * Narrowing constructor — used to create read-only views from mutable references.
         * Expansion of write permissions is forbidden at compile time.
         */
        template <ParticleField OtherWriteMask>
            requires ((WriteAccess & OtherWriteMask) == WriteAccess)
        explicit PackedParticleRef(const PackedParticleRef<ReadAccess, OtherWriteMask, Attributes, Source>& r) noexcept
            : force(r.force)
              , position(r.position)
              , velocity(r.velocity)
              , old_position(r.old_position)
              , mass(r.mass)
              , state(r.state)
              , type(r.type)
              , id(r.id)
              // , attributes(r.attributes)
        {}

        /**
          * Returns a read-only view of this SIMD block.
          */
        auto to_view() const noexcept {
            return PackedParticleRef<ReadAccess | WriteAccess, ParticleField::none, Attributes, Source>(*this);
        }

        /**
         * Loads the referenced memory block into SIMD registers (PackedParticleBuffer).
         * This is the main transition point from memory to computation.
         */
        template<MaskPolicy mask_policy = MaskPolicy::Disabled>
        PackedParticleBuffer<ReadAccess, WriteAccess, Attributes, mask_policy> load_buffer() const noexcept {
            return PackedParticleBuffer<ReadAccess, WriteAccess, Attributes, mask_policy>(*this);
        }

        // Data members with strict const-correctness
        // APRIL_NO_UNIQUE_ADDRESS ensures forbidden fields do not increase object size.
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::force> force;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::position> position;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::velocity> velocity;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::old_position> old_position;

        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::mass> mass;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::state> state;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::type> type;
        APRIL_NO_UNIQUE_ADDRESS field_t<ParticleField::id> id;

        AccessForbidden<ParticleField::attributes> attributes;
    };

    template<
        IsParticleAttributes Attributes,
        typename Source
    >
    [[nodiscard]] auto make_packed_particle_ref(Source&& source) {
        using source_t = std::remove_cvref_t<Source>;

        return PackedParticleRef<
            Source::Read,
            Source::Write,
            Attributes,
            source_t
        >(source);
    }
}
