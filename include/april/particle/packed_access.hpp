#pragma once
#include "april/base/types.hpp"
#include "april/particle/source.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"

namespace april::env {
    template<FieldMask M> struct PackedParticleView;


    template<FieldMask M, Field F, typename Source>
    constexpr auto init_packed(const Source& src) {
        if constexpr (has_field_v<M, F>) {
            return src.template get<F>();
        } else {
            return std::monostate{};
        }
    }


    //-------------------
    // PARTICLE REFERENCE
    //-------------------
    template<FieldMask M>
    struct PackedParticleRef {

        explicit PackedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            , mass        (init_packed<M, Field::mass>(source))
        {}

        PackedParticleView<M> to_view() noexcept {
            return PackedParticleView<M>(*this);
        }

        using pvec3_ref = math::Vec3Proxy<pvec3::type>;

        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::force, M> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::position, M> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::velocity, M> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::old_position, M> old_position;

        AP_NO_UNIQUE_ADDRESS field_type_t<simd::PackedRef<double>, Field::mass, M> mass;
    };



    //------------------------
    // RESTRICTED PARTICLE REF
    //------------------------
    template<FieldMask M>
    struct PackedRestrictedParticleRef {

        explicit PackedRestrictedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            , mass        (init_packed<M, Field::mass>(source))
        {}

        PackedParticleView<M> to_view() noexcept {
            return PackedParticleView<M>(*this);
        }

        using pvec3_ref = math::Vec3Proxy<pvec3::type>;
        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;

        // Force is Mutable
        AP_NO_UNIQUE_ADDRESS pvec3_ref force;

        // Others are Read-Only
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::position, M> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::velocity, M> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::old_position, M> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_d_ref, Field::mass, M> mass;
    };



    //--------------
    // PARTICLE VIEW
    //--------------
    template<FieldMask M>
    struct PackedParticleView {

        explicit PackedParticleView(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            // , mass        (init_packed<M, Field::mass>(source))
        {}

        // Copy from Mutable Ref
        template<typename RefT>
        requires (
           std::is_same_v<std::remove_cvref_t<RefT>, PackedParticleRef<M>> ||
           std::is_same_v<std::remove_cvref_t<RefT>, PackedRestrictedParticleRef<M>>
        )
        explicit PackedParticleView(const RefT& r) noexcept
           : force        (r.force)
           , position     (r.position)
           , velocity     (r.velocity)
           , old_position (r.old_position)
           // , mass         (r.mass)
        {}

        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;

        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::force, M> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::position, M> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::velocity, M> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::old_position, M> old_position;
        // AP_NO_UNIQUE_ADDRESS field_type_t<const const_d_ref, Field::mass, M> mass;
    };




}