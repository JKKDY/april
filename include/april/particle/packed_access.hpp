#pragma once
#include "april/base/types.hpp"
#include "april/particle/source.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/scalar_access.hpp"

namespace april::env {
    template<ParticleField M, IsUserData U> struct PackedParticleView;


    template<ParticleField M, ParticleField F, typename Source>
    constexpr auto init_packed(const Source& src) {
        if constexpr (has_field_v<M, F>) {
            // no de-reference here because pointers are needed to initialize packed references
            return src.template get<F>();
        } else {
            return internal::AccessForbidden<F>();
        }
    }

    //----------------
    // PACKED PARTICLE
    //----------------
    // shadow object with actual SIMD registers. Allows for manipulations without direct write backs
    template<ParticleField M>
    struct PackedParticleBuffer {
    private:
        template <typename T, ParticleField F> using field_type_t = internal::field_type_t<T, F, M>;
    public:
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3, ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3, ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3, ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3, ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<simd::Packed<double>, ParticleField::mass> mass;

        PackedParticleBuffer() = default;
        explicit PackedParticleBuffer(const auto & source) {
            if constexpr (has_field_v<M, ParticleField::position>) position = source.position;
            if constexpr (has_field_v<M, ParticleField::old_position>) old_position = source.old_position;
            if constexpr (has_field_v<M, ParticleField::velocity>) velocity = source.velocity;
            if constexpr (has_field_v<M, ParticleField::force>) force = source.force;
            if constexpr (has_field_v<M, ParticleField::mass>) mass = source.mass;
        }

        template<typename ScalarAccessor>
        requires IsScalarParticleAccessor<ScalarAccessor>
        static PackedParticleBuffer broadcast(const ScalarAccessor& scalar) {
            PackedParticleBuffer buf;

            if constexpr (has_field_v<M, ParticleField::position>) {
                buf.position.x = scalar.position.x;
                buf.position.y = scalar.position.y;
                buf.position.z = scalar.position.z;
            }
            if constexpr (has_field_v<M, ParticleField::old_position>) {
                buf.old_position.x = scalar.old_position.x;
                buf.old_position.y = scalar.old_position.y;
                buf.old_position.z = scalar.old_position.z;
            }
            if constexpr (has_field_v<M, ParticleField::velocity>) {
                buf.velocity.x = scalar.velocity.x;
                buf.velocity.y = scalar.velocity.y;
                buf.velocity.z = scalar.velocity.z;
            }
            if constexpr (has_field_v<M, ParticleField::force>) {
                buf.force.x = scalar.force.x;
                buf.force.y = scalar.force.y;
                buf.force.z = scalar.force.z;
            }
            if constexpr (has_field_v<M, ParticleField::mass>) {
                buf.mass = scalar.mass;
            }

            return buf;
        }

        template<unsigned K = 1>
        void rotate_left() {
            if constexpr (has_field_v<M, ParticleField::position>) {
                position.x = position.x.template rotate_left<K>();
                position.y = position.y.template rotate_left<K>();
                position.z = position.z.template rotate_left<K>();
            }
            if constexpr (has_field_v<M, ParticleField::old_position>) {
                old_position.x = old_position.x.template rotate_left<K>();
                old_position.y = old_position.y.template rotate_left<K>();
                old_position.z = old_position.z.template rotate_left<K>();
            }
            if constexpr (has_field_v<M, ParticleField::velocity>) {
                velocity.x = velocity.x.template rotate_left<K>();
                velocity.y = velocity.y.template rotate_left<K>();
                velocity.z = velocity.z.template rotate_left<K>();
            }
            if constexpr (has_field_v<M, ParticleField::force>) {
                force.x = force.x.template rotate_left<K>();
                force.y = force.y.template rotate_left<K>();
                force.z = force.z.template rotate_left<K>();
            }
            if constexpr (has_field_v<M, ParticleField::mass>) {
                mass = mass.template rotate_left<K>();
            }
        }

        template<unsigned K = 1>
        void rotate_right() {
            if constexpr (has_field_v<M, ParticleField::position>) {
                position.x = position.x.template rotate_right<K>();
                position.y = position.y.template rotate_right<K>();
                position.z = position.z.template rotate_right<K>();
            }
            if constexpr (has_field_v<M, ParticleField::old_position>) {
                old_position.x = old_position.x.template rotate_right<K>();
                old_position.y = old_position.y.template rotate_right<K>();
                old_position.z = old_position.z.template rotate_right<K>();
            }
            if constexpr (has_field_v<M, ParticleField::velocity>) {
                velocity.x = velocity.x.template rotate_right<K>();
                velocity.y = velocity.y.template rotate_right<K>();
                velocity.z = velocity.z.template rotate_right<K>();
            }
            if constexpr (has_field_v<M, ParticleField::force>) {
                force.x = force.x.template rotate_right<K>();
                force.y = force.y.template rotate_right<K>();
                force.z = force.z.template rotate_right<K>();
            }
            if constexpr (has_field_v<M, ParticleField::mass>) {
                mass = mass.template rotate_right<K>();
            }
        }
    };


    //-------------------
    // PARTICLE REFERENCE
    //-------------------
    template<ParticleField M, IsUserData U>
    struct PackedParticleRef {
    private:
        template <typename T, ParticleField F> using field_type_t = internal::field_type_t<T, F, M>;
        using pvec3_ref = math::Vec3Proxy<pvec3::type>;
    public:

        explicit PackedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, ParticleField::force>(source))
            , position    (init_packed<M, ParticleField::position>(source))
            , velocity    (init_packed<M, ParticleField::velocity>(source))
            , old_position(init_packed<M, ParticleField::old_position>(source))
            , mass        (init_packed<M, ParticleField::mass>(source))
        {}

        PackedParticleView<M, U> to_view() const noexcept {
            return PackedParticleView<M, U>(*this);
        }

        PackedParticleBuffer<M> load_buffer() const noexcept {
            return PackedParticleBuffer<M>(*this);
        }

        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<simd::PackedRef<double>, ParticleField::mass> mass;
    };



    //------------------------
    // RESTRICTED PARTICLE REF
    //------------------------
    template<ParticleField M, IsUserData U>
    struct PackedRestrictedParticleRef {
    private:
        template <typename T, ParticleField F> using field_type_t = internal::field_type_t<T, F, M>;
        using pvec3_ref = math::Vec3Proxy<pvec3::type>;
        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;
    public:

        explicit PackedRestrictedParticleRef(const auto& source) noexcept
            : force       (init_packed<M, ParticleField::force>(source))
            , position    (init_packed<M, ParticleField::position>(source))
            , velocity    (init_packed<M, ParticleField::velocity>(source))
            , old_position(init_packed<M, ParticleField::old_position>(source))
            , mass        (init_packed<M, ParticleField::mass>(source))
        {}

        PackedParticleView<M, U> to_view() const noexcept {
            return PackedParticleView<M, U>(*this);
        }

        PackedParticleBuffer<M> load_buffer() const noexcept {
            return PackedParticleBuffer<M>(*this);
        }

        // Force is Mutable
        AP_NO_UNIQUE_ADDRESS field_type_t<pvec3_ref, ParticleField::force> force;

        // Others are Read-Only
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_d_ref, ParticleField::mass> mass;
    };



    //--------------
    // PARTICLE VIEW
    //--------------
    template<ParticleField M, IsUserData U>
    struct PackedParticleView {
    private:
        template <typename T, ParticleField F> using field_type_t = internal::field_type_t<T, F, M>;
        using const_pvec3_ref = math::Vec3Proxy<const pvec3::type>;
        using const_d_ref = simd::PackedRef<const double>;
    public:

        explicit PackedParticleView(const auto& source) noexcept
            : force       (init_packed<M, ParticleField::force>(source))
            , position    (init_packed<M, ParticleField::position>(source))
            , velocity    (init_packed<M, ParticleField::velocity>(source))
            , old_position(init_packed<M, ParticleField::old_position>(source))
            , mass        (init_packed<M, ParticleField::mass>(source))
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

        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_pvec3_ref, ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS field_type_t<const const_d_ref, ParticleField::mass> mass;
    };


    //---------
    // CONCEPTS
    //---------
    template<typename T> struct is_packed_buffer_impl         : std::false_type {};
    template<typename T> struct is_packed_ref_impl            : std::false_type {};
    template<typename T> struct is_packed_restricted_ref_impl : std::false_type {};
    template<typename T> struct is_packed_view_impl           : std::false_type {};

    // Specialization for Accessors
    template<auto M>             struct is_packed_buffer_impl<PackedParticleBuffer<M>> : std::true_type {};
    template<auto M, typename U> struct is_packed_ref_impl<PackedParticleRef<M, U>> : std::true_type {};
    template<auto M, typename U> struct is_packed_restricted_ref_impl<PackedRestrictedParticleRef<M, U>> : std::true_type {};
    template<auto M, typename U> struct is_packed_view_impl<PackedParticleView<M, U>> : std::true_type {};

	// Concepts
    template<typename T>
    concept IsPackedParticleBuffer = is_packed_buffer_impl<std::remove_cvref_t<T>>::value;

    template<typename T>
    concept IsPackedParticleRef = is_packed_ref_impl<std::remove_cvref_t<T>>::value;

    template<typename T>
    concept IsPackedRestrictedParticleRef = is_packed_restricted_ref_impl<std::remove_cvref_t<T>>::value;

    template<typename T>
   concept IsPackedParticleView = is_packed_view_impl<std::remove_cvref_t<T>>::value;

    template<typename T>
    concept IsPackedParticleAccessor =
        IsPackedParticleBuffer<T> ||
        IsPackedParticleRef<T> ||
        IsPackedRestrictedParticleRef<T> ||
        IsPackedParticleView<T>;

    template<typename T>
    concept IsAnyParticleAccessor = IsScalarParticleAccessor<T> || IsPackedParticleAccessor<T>;
}




