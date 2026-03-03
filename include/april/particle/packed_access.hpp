#pragma once
#include "april/base/types.hpp"
#include "april/particle/source.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/scalar_access.hpp"

namespace april::particle::internal {
    // fwd declaration
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedBufferView;
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes U>
    struct PackedParticleRef;


    //-----------------------
    // PACKED PARTICLE BUFFER
    //-----------------------
    // if in ReadMask or WriteMask return type T else return Poison


    // shadow object with actual SIMD registers. Allows for manipulations without direct write backs
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedParticleBuffer {
    private:
        template <typename T, ParticleField F>
          using buffer_field_t = std::conditional_t<
              has_field_v<ReadMask | WriteMask, F>,
              T,
              AccessForbidden<F>
          >;

        template <ParticleField F>
        using pvec3_t = buffer_field_t<pvec3, F>;

        template <ParticleField F>
        using scalar_t = buffer_field_t<simd::Packed<double>, F>;
    public:
        static constexpr ParticleField RWMask = ReadMask & WriteMask;  // read & write
        static constexpr ParticleField WOMask = WriteMask & ~ReadMask; // write only
        static constexpr ParticleField ROMask = ReadMask & ~WriteMask; // read only

        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS scalar_t<ParticleField::mass> mass;

        PackedParticleBuffer() = default;

        // if in read mask read from memory else if (only) in write mask zero initialize
        template <typename attr>
        explicit PackedParticleBuffer(const PackedParticleRef<ReadMask, WriteMask, attr>& source) {
            if constexpr (has_field_v<ReadMask, ParticleField::position>) position = source.position;
            if constexpr (has_field_v<ReadMask, ParticleField::old_position>) old_position = source.old_position;
            if constexpr (has_field_v<ReadMask, ParticleField::velocity>) velocity = source.velocity;
            if constexpr (has_field_v<ReadMask, ParticleField::force>) force = source.force;
            if constexpr (has_field_v<ReadMask, ParticleField::mass>) mass = source.mass;

            // Write-Only fields: Zero-initialize for pure delta accumulation (necessary for symmetric batches)
            if constexpr (has_field_v<WOMask, ParticleField::position>) position = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) old_position = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) velocity = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::force>) force = pvec3(0.0);
            if constexpr (has_field_v<WOMask, ParticleField::mass>) mass = 0.0;
        }

        template <typename ScalarAccessor>
            requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        explicit PackedParticleBuffer(const ScalarAccessor& scalar) {
            if constexpr (has_field_v<ReadMask, ParticleField::position>) {
                position.x = scalar.position.x;
                position.y = scalar.position.y;
                position.z = scalar.position.z;
            }
            else if constexpr (has_field_v<WOMask, ParticleField::position>) {
                position = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::old_position>) {
                old_position.x = scalar.old_position.x;
                old_position.y = scalar.old_position.y;
                old_position.z = scalar.old_position.z;
            }
            else if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                old_position = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::velocity>) {
                velocity.x = scalar.velocity.x;
                velocity.y = scalar.velocity.y;
                velocity.z = scalar.velocity.z;
            }
            else if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                velocity = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::force>) {
                force.x = scalar.force.x;
                force.y = scalar.force.y;
                force.z = scalar.force.z;
            }
            else if constexpr (has_field_v<WOMask, ParticleField::force>) {
                force = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::mass>) {
                mass = scalar.mass;
            }
            else if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass = 0.0;
            }
        }

        // export as view
        AP_FORCE_INLINE PackedBufferView<ReadMask, WriteMask> to_view() {
            return PackedBufferView(*this);
        }

        // Trivial in-place rotations
        template <unsigned K = 1>
        AP_FORCE_INLINE void rotate_left() {
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::position>) {
                position.x = position.x.template rotate_left<K>();
                position.y = position.y.template rotate_left<K>();
                position.z = position.z.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::old_position>) {
                old_position.x = old_position.x.template rotate_left<K>();
                old_position.y = old_position.y.template rotate_left<K>();
                old_position.z = old_position.z.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::velocity>) {
                velocity.x = velocity.x.template rotate_left<K>();
                velocity.y = velocity.y.template rotate_left<K>();
                velocity.z = velocity.z.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::force>) {
                force.x = force.x.template rotate_left<K>();
                force.y = force.y.template rotate_left<K>();
                force.z = force.z.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::mass>) {
                mass = mass.template rotate_left<K>();
            }
        }

        template <unsigned K = 1>
        AP_FORCE_INLINE void rotate_right() {
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::position>) {
                position.x = position.x.template rotate_right<K>();
                position.y = position.y.template rotate_right<K>();
                position.z = position.z.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::old_position>) {
                old_position.x = old_position.x.template rotate_right<K>();
                old_position.y = old_position.y.template rotate_right<K>();
                old_position.z = old_position.z.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::velocity>) {
                velocity.x = velocity.x.template rotate_right<K>();
                velocity.y = velocity.y.template rotate_right<K>();
                velocity.z = velocity.z.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::force>) {
                force.x = force.x.template rotate_right<K>();
                force.y = force.y.template rotate_right<K>();
                force.z = force.z.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::mass>) {
                mass = mass.template rotate_right<K>();
            }
        }


        // Unmasked SIMD Write-Back (For full chunks)
        template <typename Attr>
        AP_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr>& packed_ref) const {
            // Write-Only fields use additive accumulation (preserves base state)

            // POSITION
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                packed_ref.position += position;
            } else if constexpr (has_field_v<RWMask, ParticleField::position>) {
                packed_ref.position = position;
            }

            // OLD POSITION
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                packed_ref.old_position += old_position;
            } else if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                packed_ref.old_position = old_position;
            }

            // VELOCITY
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                packed_ref.velocity += velocity;
            } else if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                packed_ref.velocity = velocity;
            }

            // FORCE
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                packed_ref.force += force;
            } else if constexpr (has_field_v<RWMask, ParticleField::force>) {
                packed_ref.force = force;
            }

            // MASS
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                packed_ref.mass += mass;
            } else if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                packed_ref.mass = mass;
            }
        }

        // Masked SIMD Write-Back (For Tail Chunks)
        template <typename Attr, typename MaskT>
        AP_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr>& packed_ref, const MaskT & mask) const {
            const packed null = 0.0;

            // POSITION
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                packed_ref.position.x += select(mask, position.x, null);
                packed_ref.position.y += select(mask, position.y, null);
                packed_ref.position.z += select(mask, position.z, null);
            } else if constexpr (has_field_v<RWMask, ParticleField::position>) {
                packed_ref.position.x = select(mask, position.x, packed_ref.position.x);
                packed_ref.position.y = select(mask, position.y, packed_ref.position.y);
                packed_ref.position.z = select(mask, position.z, packed_ref.position.z);
            }

            // OLD POSITION
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                packed_ref.old_position.x += select(mask, old_position.x, null);
                packed_ref.old_position.y += select(mask, old_position.y, null);
                packed_ref.old_position.z += select(mask, old_position.z, null);
            } else if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                packed_ref.old_position.x = select(mask, old_position.x, packed_ref.old_position.x);
                packed_ref.old_position.y = select(mask, old_position.y, packed_ref.old_position.y);
                packed_ref.old_position.z = select(mask, old_position.z, packed_ref.old_position.z);
            }

            // VELOCITY
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                packed_ref.velocity.x += select(mask, velocity.x, null);
                packed_ref.velocity.y += select(mask, velocity.y, null);
                packed_ref.velocity.z += select(mask, velocity.z, null);
            } else if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                packed_ref.velocity.x = select(mask, velocity.x, packed_ref.velocity.x);
                packed_ref.velocity.y = select(mask, velocity.y, packed_ref.velocity.y);
                packed_ref.velocity.z = select(mask, velocity.z, packed_ref.velocity.z);
            }

            // FORCE
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                packed_ref.force.x += select(mask, force.x, null);
                packed_ref.force.y += select(mask, force.y, null);
                packed_ref.force.z += select(mask, force.z, null);
            } else if constexpr (has_field_v<RWMask, ParticleField::force>) {
                packed_ref.force.x = select(mask, force.x, packed_ref.force.x);
                packed_ref.force.y = select(mask, force.y, packed_ref.force.y);
                packed_ref.force.z = select(mask, force.z, packed_ref.force.z);
            }

            // MASS
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                packed_ref.mass += select(mask, mass, null);
            } else if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                packed_ref.mass = select(mask, mass, packed_ref.mass);
            }
        }

        // Unmasked Scalar Reduction (For Broadcast Buffers)
        template <typename ScalarAccessor>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p) const {

            // POSITION
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                p.position.x += position.x.reduce_add();
                p.position.y += position.y.reduce_add();
                p.position.z += position.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::position>) {
                static_assert(false,
                    "FATAL: Cannot reduce a Read-Write 'position' field from a SIMD register to a scalar. "
                    "This implies divergent state updates across lanes, which is mathematically invalid.");
            }

            // OLD POSITION
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                p.old_position.x += old_position.x.reduce_add();
                p.old_position.y += old_position.y.reduce_add();
                p.old_position.z += old_position.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                static_assert(false,
                    "FATAL: Cannot reduce a Read-Write 'old_position' field from a SIMD register to a scalar.");
            }

            // VELOCITY
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                p.velocity.x += velocity.x.reduce_add();
                p.velocity.y += velocity.y.reduce_add();
                p.velocity.z += velocity.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                static_assert(false,
                    "FATAL: Cannot reduce a Read-Write 'velocity' field from a SIMD register to a scalar.");
            }

            // FORCE
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                p.force.x += force.x.reduce_add();
                p.force.y += force.y.reduce_add();
                p.force.z += force.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::force>) {
                static_assert(false,
                    "FATAL: Cannot reduce a Read-Write 'force' field from a SIMD register to a scalar.");
            }

            // MASS
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                p.mass += mass.reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                static_assert(false,
                    "FATAL: Cannot reduce a Read-Write 'mass' field from a SIMD register to a scalar.");
            }
        }

        // Masked Scalar Reduction (For partial tail logic vs broadcasted scalars)
        template <typename ScalarAccessor, typename MaskT>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p, const MaskT& mask) const {
            const packed null = 0.0;

            // POSITION
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                p.position.x += select(mask, position.x, null).reduce_add();
                p.position.y += select(mask, position.y, null).reduce_add();
                p.position.z += select(mask, position.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::position>) {
                static_assert(false,
                    "FATAL: Cannot perform masked reduction on a Read-Write 'position' field.");
            }

            // OLD POSITION
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                p.old_position.x += select(mask, old_position.x, null).reduce_add();
                p.old_position.y += select(mask, old_position.y, null).reduce_add();
                p.old_position.z += select(mask, old_position.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                static_assert(false,
                    "FATAL: Cannot perform masked reduction on a Read-Write 'old_position' field.");
            }

            // VELOCITY
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                p.velocity.x += select(mask, velocity.x, null).reduce_add();
                p.velocity.y += select(mask, velocity.y, null).reduce_add();
                p.velocity.z += select(mask, velocity.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                static_assert(false,
                    "FATAL: Cannot perform masked reduction on a Read-Write 'velocity' field.");
            }

            // FORCE
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                p.force.x += select(mask, force.x, null).reduce_add();
                p.force.y += select(mask, force.y, null).reduce_add();
                p.force.z += select(mask, force.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::force>) {
                static_assert(false,
                    "FATAL: Cannot perform masked reduction on a Read-Write 'force' field.");
            }

            // MASS
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                p.mass += select(mask, mass, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                static_assert(false,
                    "FATAL: Cannot perform masked reduction on a Read-Write 'mass' field.");
            }
        }
    };


    //------------
    // BUFFER VIEW
    //------------
    template <ParticleField F, ParticleField WriteMask, typename T>
    using view_ref_t = std::conditional_t<
        std::is_same_v<T, AccessForbidden<F>>, // if it's poison, keep it as poison (by value)
        T,
        std::conditional_t<has_field_v<WriteMask, F>, T&, const T&>
        // if it's valid and in WriteMask, make it a mutable ref else a const ref
    >;


    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedBufferView {
        using Buffer = PackedParticleBuffer<ReadMask, WriteMask>;

    private:
        template <ParticleField F, typename T>
        using ref_t = view_ref_t<F, WriteMask, T>;

    public:
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::position, decltype(Buffer::position)> position;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::old_position, decltype(Buffer::old_position)> old_position;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::velocity, decltype(Buffer::velocity)> velocity;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::force, decltype(Buffer::force)> force;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::mass, decltype(Buffer::mass)> mass;

        // Directly maps buffer fields to references or poison copies. Zero branching.
        AP_FORCE_INLINE explicit PackedBufferView(Buffer& buf)
            : position(buf.position),
              old_position(buf.old_position),
              velocity(buf.velocity),
              force(buf.force),
              mass(buf.mass) {}
    };



    //--------------------
    // PACKED PARTICLE REF
    //--------------------
    // helper for switching between data pointers and poison
    template <ParticleField ReadMask, ParticleField WriteMask, ParticleField F, typename Source>
    constexpr auto init_packed(const Source& src) {
        if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, F>) {
            // no de-reference here because pointers are needed to initialize packed references
            return src.template get<F>();
        }
        else {
            return AccessForbidden<F>();
        }
    }

    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes U>
    struct PackedParticleRef {
    private:
        template <typename MutT, typename ConstT, ParticleField F>
        using field_t = field_access_t<MutT, ConstT, F, ReadMask, WriteMask>;

        using MutVec3Ref = math::Vec3Proxy<pvec3::type>;
        using ConstVec3Ref = math::Vec3Proxy<const pvec3::type>;

    public:
        explicit PackedParticleRef(const auto& source) noexcept
            : force(internal::init_packed<ReadMask, WriteMask, ParticleField::force>(source))
              , position(internal::init_packed<ReadMask, WriteMask, ParticleField::position>(source))
              , velocity(internal::init_packed<ReadMask, WriteMask, ParticleField::velocity>(source))
              , old_position(internal::init_packed<ReadMask, WriteMask, ParticleField::old_position>(source))
              , mass(internal::init_packed<ReadMask, WriteMask, ParticleField::mass>(source)) {}

        // Cross-constructor for narrowing write permissions (e.g. creating a view)
        template <ParticleField OtherWriteMask>
            requires ((WriteMask & OtherWriteMask) == WriteMask)
        explicit PackedParticleRef(const PackedParticleRef<ReadMask, OtherWriteMask, U>& r) noexcept
            : force(r.force)
              , position(r.position)
              , velocity(r.velocity)
              , old_position(r.old_position)
              , mass(r.mass) {}

        // A View is just a PackedParticleRef with no write permissions
        auto to_view() const noexcept {
            return PackedParticleRef<ReadMask, ParticleField::none, U>(*this);
        }

        // Spits out the mutable buffer we just designed
        PackedParticleBuffer<ReadMask, WriteMask> load_buffer() const noexcept {
            return PackedParticleBuffer<ReadMask, WriteMask>(*this);
        }

        // Write the buffer back to memory, strictly respecting the WriteMask
        AP_FORCE_INLINE void update(const PackedParticleBuffer<ReadMask, WriteMask>& buf) noexcept {
            if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::position>) {
                position = buf.position;
            }
            if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::old_position>) {
                old_position = buf.old_position;
            }
            if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::velocity>) {
                velocity = buf.velocity;
            }
            if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::force>) {
                force = buf.force;
            }
            if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::mass>) {
                mass = buf.mass;
            }
        }

        // Data members with strict const-correctness
        AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref, ConstVec3Ref, ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref, ConstVec3Ref, ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref, ConstVec3Ref, ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref, ConstVec3Ref, ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS field_t<simd::PackedRef<double>, simd::PackedRef<const double>, ParticleField::mass> mass;
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
    template <ParticleField RM, ParticleField WM>
    struct is_packed_buffer_impl<PackedParticleBuffer<RM, WM>> : std::true_type {};

    template <ParticleField RM, ParticleField WM, typename U>
    struct is_packed_ref_impl<PackedParticleRef<RM, WM, U>> : std::true_type {};

    template <ParticleField RM, ParticleField WM>
    struct is_buffer_view_impl<PackedBufferView<RM, WM>> : std::true_type {};
} // namespace april::particle::internal


namespace april::particle {
    //---------
    // CONCEPTS
    //---------
    // TODO review which types to expose
    template <typename T>
    concept IsPackedParticleBuffer = internal::is_packed_buffer_impl<std::remove_cvref_t<T>>::value;

    template <typename T>
    concept IsPackedParticleRef = internal::is_packed_ref_impl<std::remove_cvref_t<T>>::value;

    template <typename T>
    concept IsPackedParticleView = internal::is_buffer_view_impl<std::remove_cvref_t<T>>::value;

    template <typename T>
    concept IsPackedParticleAccessor =
        IsPackedParticleBuffer<T> ||
        IsPackedParticleRef<T> ||
        IsPackedParticleView<T>;

    template <typename T>
    concept IsAnyParticleAccessor = IsScalarParticleAccessor<T> || IsPackedParticleAccessor<T>;
}
