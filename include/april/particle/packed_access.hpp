#pragma once
#include "april/base/types.hpp"
#include "april/particle/source.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"

namespace april::env {
    template<FieldMask M, IsUserData U> struct PackedParticleView;


    template<FieldMask M, Field F, typename Source>
    constexpr auto init_packed(const Source& src) {
        if constexpr (has_field_v<M, F>) {
            return src.template get<F>();
        } else {
            return std::monostate{};
        }
    }


    //----------------
    // PACKED PARTICLE
    //----------------
    // shadow object with actual SIMD registers. Allows for manipulations without direct write backs
    template<FieldMask M>
    struct PackedParticleBuffer {
        field_type_t<pvec3, Field::position, M> position;
        field_type_t<pvec3, Field::old_position, M> old_position;
        field_type_t<pvec3, Field::velocity, M> velocity;
        field_type_t<pvec3, Field::force, M> force;
        field_type_t<simd::Packed<double>, Field::mass, M> mass;

        explicit PackedParticleBuffer(const auto & source) {
            if constexpr (has_field_v<M, Field::position>) position = source.position;
            if constexpr (has_field_v<M, Field::old_position>) old_position = source.old_position;
            if constexpr (has_field_v<M, Field::velocity>) velocity = source.velocity;
            if constexpr (has_field_v<M, Field::force>) force = source.force;
            if constexpr (has_field_v<M, Field::mass>) mass = source.mass;
        }

        void rotate_left() {
            if constexpr (has_field_v<M, Field::position>) {
                position.x = position.x.rotate_left();
                position.y = position.y.rotate_left();
                position.z = position.z.rotate_left();
            }
            if constexpr (has_field_v<M, Field::old_position>) {
                old_position.x = old_position.x.rotate_left();
                old_position.y = old_position.y.rotate_left();
                old_position.z = old_position.z.rotate_left();
            }
            if constexpr (has_field_v<M, Field::velocity>) {
                velocity.x = velocity.x.rotate_left();
                velocity.y = velocity.y.rotate_left();
                velocity.z = velocity.z.rotate_left();
            }
            if constexpr (has_field_v<M, Field::force>) {
                force.x = force.x.rotate_left();
                force.y = force.y.rotate_left();
                force.z = force.z.rotate_left();
            }
            if constexpr (has_field_v<M, Field::mass>) mass = mass.rotate_left();
        }

        void rotate_right() {
            if constexpr (has_field_v<M, Field::position>) {
                position.x = position.x.rotate_right();
                position.y = position.y.rotate_right();
                position.z = position.z.rotate_right();
            }
            if constexpr (has_field_v<M, Field::old_position>) {
                old_position.x = old_position.x.rotate_right();
                old_position.y = old_position.y.rotate_right();
                old_position.z = old_position.z.rotate_right();
            }
            if constexpr (has_field_v<M, Field::velocity>) {
                velocity.x = velocity.x.rotate_right();
                velocity.y = velocity.y.rotate_right();
                velocity.z = velocity.z.rotate_right();
            }
            if constexpr (has_field_v<M, Field::force>) {
                force.x = force.x.rotate_right();
                force.y = force.y.rotate_right();
                force.z = force.z.rotate_right();
            }
            if constexpr (has_field_v<M, Field::mass>) mass = mass.rotate_right();
        }
    };


    //-------------------
    // PARTICLE REFERENCE
    //-------------------
    template<FieldMask M, IsUserData U>
    struct PackedParticleRef {
        using pvec3_ref = math::Vec3Proxy<pvec3::type>;

        explicit PackedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            , mass        (init_packed<M, Field::mass>(source))
        {}

        PackedParticleView<M, U> to_view() const noexcept {
            return PackedParticleView<M, U>(*this);
        }

        PackedParticleBuffer<M> load_buffer() const noexcept {
            return PackedParticleBuffer<M>(*this);
        }


        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::force, M> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::position, M> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::velocity, M> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, Field::old_position, M> old_position;

        AP_NO_UNIQUE_ADDRESS field_type_t<simd::PackedRef<double>, Field::mass, M> mass;
    };



    //------------------------
    // RESTRICTED PARTICLE REF
    //------------------------
    template<FieldMask M, IsUserData U>
    struct PackedRestrictedParticleRef {
        using pvec3_ref = math::Vec3Proxy<pvec3::type>;
        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;

        explicit PackedRestrictedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            , mass        (init_packed<M, Field::mass>(source))
        {}

        PackedParticleView<M, U> to_view() const noexcept {
            return PackedParticleView<M, U>(*this);
        }

        PackedParticleBuffer<M> load_buffer() const noexcept {
            return PackedParticleBuffer<M>(*this);
        }

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
    template<FieldMask M, IsUserData U>
    struct PackedParticleView {
        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;

        explicit PackedParticleView(const auto& source) noexcept
            : force       (init_packed<M, Field::force>(source))
            , position    (init_packed<M, Field::position>(source))
            , velocity    (init_packed<M, Field::velocity>(source))
            , old_position(init_packed<M, Field::old_position>(source))
            , mass        (init_packed<M, Field::mass>(source))
        {}

        // Copy from Mutable Ref
        template<typename RefT>
        requires (
           std::is_same_v<std::remove_cvref_t<RefT>, PackedParticleRef<M, U>> ||
           std::is_same_v<std::remove_cvref_t<RefT>, PackedRestrictedParticleRef<M, U>>
        )
        explicit PackedParticleView(const RefT& r) noexcept
           : force        (r.force)
           , position     (r.position)
           , velocity     (r.velocity)
           , old_position (r.old_position)
           , mass         (r.mass)
        {}

        PackedParticleBuffer<M> load_buffer() const noexcept {
            return PackedParticleBuffer<M>(*this);
        }

        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::force, M> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::position, M> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::velocity, M> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, Field::old_position, M> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_d_ref, Field::mass, M> mass;
    };


    //


}