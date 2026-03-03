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


    //-----------------------
    // PACKED PARTICLE BUFFER
    //-----------------------
    // if in ReadMask or WriteMask return type T else return Poison
    template <typename T, ParticleField F, ParticleField ReadMask, ParticleField WriteMask>
    using buffer_field_t = std::conditional_t<
        has_field_v<ReadMask | WriteMask, F>,
        T,
        AccessForbidden<F>
    >;

    // shadow object with actual SIMD registers. Allows for manipulations without direct write backs
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedParticleBuffer {
    private:
        template <ParticleField F>
        using pvec3_t = buffer_field_t<pvec3, F, ReadMask, WriteMask>;

        template <ParticleField F>
        using scalar_t = buffer_field_t<simd::Packed<double>, F, ReadMask, WriteMask>;

    public:
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS scalar_t<ParticleField::mass> mass;

        PackedParticleBuffer() = default;

        // if in read mask read from memory else if (only) in write mask zero initialize
        explicit PackedParticleBuffer(const auto& source) {
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::position>) position =
                source.position;
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::old_position>)
                old_position = source.old_position;
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::velocity>) velocity =
                source.velocity;
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::force>) force = source.
                force;
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::mass>) mass = source.
                mass;
        }

        // Broadcast from a scalar particle
        template <typename ScalarAccessor>
            requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        static PackedParticleBuffer broadcast(const ScalarAccessor& scalar) {
            PackedParticleBuffer buf;

            constexpr auto Mask = ReadMask | WriteMask;

            if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>) {
                buf.position.x = scalar.position.x;
                buf.position.y = scalar.position.y;
                buf.position.z = scalar.position.z;
            }
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>) {
                buf.old_position.x = scalar.old_position.x;
                buf.old_position.y = scalar.old_position.y;
                buf.old_position.z = scalar.old_position.z;
            }
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>) {
                buf.velocity.x = scalar.velocity.x;
                buf.velocity.y = scalar.velocity.y;
                buf.velocity.z = scalar.velocity.z;
            }
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>) {
                buf.force.x = scalar.force.x;
                buf.force.y = scalar.force.y;
                buf.force.z = scalar.force.z;
            }
            if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>) {
                buf.mass = scalar.mass;
            }

            return buf;
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
    };


    //-------------------
    // PACKED BUFFER VIEW
    //-------------------
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


    //-------------------
    // PACKED ACCUMULATOR
    //-------------------
    template <ParticleField RWMask>
    struct DeltaAccumulator {
        using BufferT = PackedParticleBuffer<RWMask, RWMask>;

        BufferT pristine; // holds original data
        BufferT acc; // holds accumulated difference to original data

        AP_FORCE_INLINE DeltaAccumulator() {
            reset();
        }

        template <ParticleField RM, ParticleField WM>
        AP_FORCE_INLINE void save_pristine(const PackedParticleBuffer<RM, WM>& buf) {
            if constexpr (has_field_v<RWMask, ParticleField::position>) pristine.position = buf.position;
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) pristine.old_position = buf.old_position;
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) pristine.velocity = buf.velocity;
            if constexpr (has_field_v<RWMask, ParticleField::force>) pristine.force = buf.force;
            if constexpr (has_field_v<RWMask, ParticleField::mass>) pristine.mass = buf.mass;
        }

        void reset() {
            if constexpr (has_field_v<RWMask, ParticleField::position>) acc.position = {0, 0, 0};
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) acc.old_position = {0, 0, 0};
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) acc.velocity = {0, 0, 0};
            if constexpr (has_field_v<RWMask, ParticleField::force>) acc.force = {0, 0, 0};
            if constexpr (has_field_v<RWMask, ParticleField::mass>) acc.mass = 0.0;
        }

        template <ParticleField RM, ParticleField WM>
        AP_FORCE_INLINE void extract_and_restore(PackedParticleBuffer<RM, WM>& buf) {
            if constexpr (has_field_v<RWMask, ParticleField::position>) {
                acc.position += (buf.position - pristine.position);
                buf.position = pristine.position;
            }
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                acc.old_position += (buf.old_position - pristine.old_position);
                buf.old_position = pristine.old_position;
            }
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                acc.velocity += (buf.velocity - pristine.velocity);
                buf.velocity = pristine.velocity;
            }
            if constexpr (has_field_v<RWMask, ParticleField::force>) {
                acc.force += (buf.force - pristine.force);
                buf.force = pristine.force;
            }
            if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                acc.mass += (buf.mass - pristine.mass);
                buf.mass = pristine.mass;
            }
        }

        template <unsigned K = 1>
        AP_FORCE_INLINE void rotate_right() {
            acc.template rotate_right<K>();
        }

        template <unsigned K = 1>
        AP_FORCE_INLINE void rotate_left() {
            acc.template rotate_left<K>();
        }

        template <ParticleField RM, ParticleField WM>
        AP_FORCE_INLINE void apply_to(PackedParticleBuffer<RM, WM>& buf) const {
            if constexpr (has_field_v<RWMask, ParticleField::position>) buf.position += acc.position;
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) buf.old_position += acc.old_position;
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) buf.velocity += acc.velocity;
            if constexpr (has_field_v<RWMask, ParticleField::force>) buf.force += acc.force;
            if constexpr (has_field_v<RWMask, ParticleField::mass>) buf.mass += acc.mass;
        }
    };


    //------------------
    // BROADCAST REDUCER
    //------------------
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct BroadcastReducer {
        static constexpr ParticleField RWMask = ReadMask & WriteMask;
        static constexpr ParticleField WOMask = WriteMask & ~ReadMask; // Write-Only Mask

        PackedParticleBuffer<ReadMask, WriteMask> buffer;

        // We only need to store base values if a field is explicitly Read+Write
        AP_NO_UNIQUE_ADDRESS PackedParticleBuffer<RWMask, ParticleField::none> pristine;

        template <typename ScalarAccessor>
            requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        AP_FORCE_INLINE explicit BroadcastReducer(const ScalarAccessor& scalar) {
            // Read field: load data into buffer
            // Write only field: load a zero vector for accumulation
            // Read + Write field: store pristine data

            // POSITION
            if constexpr (has_field_v<ReadMask, ParticleField::position>)
                buffer.position = pvec3(scalar.position.x, scalar.position.y, scalar.position.z);
            else if constexpr (has_field_v<WOMask, ParticleField::position>)
                buffer.position = pvec3(0.0);

            // OLD POSITION
            if constexpr (has_field_v<ReadMask, ParticleField::old_position>)
                buffer.old_position = pvec3(scalar.old_position.x, scalar.old_position.y, scalar.old_position.z);
            else if constexpr (has_field_v<WOMask, ParticleField::old_position>)
                buffer.old_position = pvec3(0.0);

            // VELOCITY
            if constexpr (has_field_v<ReadMask, ParticleField::velocity>)
                buffer.velocity = pvec3(scalar.velocity.x, scalar.velocity.y, scalar.velocity.z);
            else if constexpr (has_field_v<WOMask, ParticleField::velocity>)
                buffer.velocity = pvec3(0.0);

            // FORCE
            if constexpr (has_field_v<ReadMask, ParticleField::force>)
                buffer.force = pvec3(scalar.force.x, scalar.force.y, scalar.force.z);
            else if constexpr (has_field_v<WOMask, ParticleField::force>)
                buffer.force = pvec3(0.0);

            // MASS
            if constexpr (has_field_v<ReadMask, ParticleField::mass>)
                buffer.mass = scalar.mass;
            else if constexpr (has_field_v<WOMask, ParticleField::mass>)
                buffer.mass = 0.0;

            // Save the base state ONLY for fields we are reading and writing
            if constexpr (has_field_v<RWMask, ParticleField::position>) pristine.position = buffer.position;
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) pristine.old_position = buffer.old_position;
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) pristine.velocity = buffer.velocity;
            if constexpr (has_field_v<RWMask, ParticleField::force>) pristine.force = buffer.force;
            if constexpr (has_field_v<RWMask, ParticleField::mass>) pristine.mass = buffer.mass;
        }

        AP_FORCE_INLINE auto to_view() {
            return buffer.to_view();
        }

        // Unmasked Reduction
        template <typename ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p) const {
            if constexpr (has_field_v<WriteMask, ParticleField::position>) {
                if constexpr (has_field_v<RWMask, ParticleField::position>) {
                    p.position.x += (buffer.position.x - pristine.position.x).reduce_add();
                    p.position.y += (buffer.position.y - pristine.position.y).reduce_add();
                    p.position.z += (buffer.position.z - pristine.position.z).reduce_add();
                }
                else {
                    p.position.x += buffer.position.x.reduce_add();
                    p.position.y += buffer.position.y.reduce_add();
                    p.position.z += buffer.position.z.reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::old_position>) {
                if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                    p.old_position.x += (buffer.old_position.x - pristine.old_position.x).reduce_add();
                    p.old_position.y += (buffer.old_position.y - pristine.old_position.y).reduce_add();
                    p.old_position.z += (buffer.old_position.z - pristine.old_position.z).reduce_add();
                }
                else {
                    p.old_position.x += buffer.old_position.x.reduce_add();
                    p.old_position.y += buffer.old_position.y.reduce_add();
                    p.old_position.z += buffer.old_position.z.reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::velocity>) {
                if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                    p.velocity.x += (buffer.velocity.x - pristine.velocity.x).reduce_add();
                    p.velocity.y += (buffer.velocity.y - pristine.velocity.y).reduce_add();
                    p.velocity.z += (buffer.velocity.z - pristine.velocity.z).reduce_add();
                }
                else {
                    p.velocity.x += buffer.velocity.x.reduce_add();
                    p.velocity.y += buffer.velocity.y.reduce_add();
                    p.velocity.z += buffer.velocity.z.reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::force>) {
                if constexpr (has_field_v<RWMask, ParticleField::force>) {
                    p.force.x += (buffer.force.x - pristine.force.x).reduce_add();
                    p.force.y += (buffer.force.y - pristine.force.y).reduce_add();
                    p.force.z += (buffer.force.z - pristine.force.z).reduce_add();
                }
                else {
                    p.force.x += buffer.force.x.reduce_add();
                    p.force.y += buffer.force.y.reduce_add();
                    p.force.z += buffer.force.z.reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::mass>) {
                if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                    p.mass += (buffer.mass - pristine.mass).reduce_add();
                }
                else {
                    p.mass += buffer.mass.reduce_add();
                }
            }
        }

        // Masked Reduction (For partial tail logic)
        // Masked Reduction (For partial tail logic)
        template <typename ScalarAccessor, typename MaskT>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p, const MaskT& mask) const {
            const packed null = 0.0;

            if constexpr (has_field_v<WriteMask, ParticleField::position>) {
                if constexpr (has_field_v<RWMask, ParticleField::position>) {
                    p.position.x += select(mask, buffer.position.x - pristine.position.x, null).reduce_add();
                    p.position.y += select(mask, buffer.position.y - pristine.position.y, null).reduce_add();
                    p.position.z += select(mask, buffer.position.z - pristine.position.z, null).reduce_add();
                }
                else {
                    p.position.x += select(mask, buffer.position.x, null).reduce_add();
                    p.position.y += select(mask, buffer.position.y, null).reduce_add();
                    p.position.z += select(mask, buffer.position.z, null).reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::old_position>) {
                if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                    p.old_position.x += select(mask, buffer.old_position.x - pristine.old_position.x, null).reduce_add();
                    p.old_position.y += select(mask, buffer.old_position.y - pristine.old_position.y, null).reduce_add();
                    p.old_position.z += select(mask, buffer.old_position.z - pristine.old_position.z, null).reduce_add();
                }
                else {
                    p.old_position.x += select(mask, buffer.old_position.x, null).reduce_add();
                    p.old_position.y += select(mask, buffer.old_position.y, null).reduce_add();
                    p.old_position.z += select(mask, buffer.old_position.z, null).reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::velocity>) {
                if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                    p.velocity.x += select(mask, buffer.velocity.x - pristine.velocity.x, null).reduce_add();
                    p.velocity.y += select(mask, buffer.velocity.y - pristine.velocity.y, null).reduce_add();
                    p.velocity.z += select(mask, buffer.velocity.z - pristine.velocity.z, null).reduce_add();
                }
                else {
                    p.velocity.x += select(mask, buffer.velocity.x, null).reduce_add();
                    p.velocity.y += select(mask, buffer.velocity.y, null).reduce_add();
                    p.velocity.z += select(mask, buffer.velocity.z, null).reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::force>) {
                if constexpr (has_field_v<RWMask, ParticleField::force>) {
                    p.force.x += select(mask, buffer.force.x - pristine.force.x, null).reduce_add();
                    p.force.y += select(mask, buffer.force.y - pristine.force.y, null).reduce_add();
                    p.force.z += select(mask, buffer.force.z - pristine.force.z, null).reduce_add();
                }
                else {
                    p.force.x += select(mask, buffer.force.x, null).reduce_add();
                    p.force.y += select(mask, buffer.force.y, null).reduce_add();
                    p.force.z += select(mask, buffer.force.z, null).reduce_add();
                }
            }
            if constexpr (has_field_v<WriteMask, ParticleField::mass>) {
                if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                    p.mass += select(mask, buffer.mass - pristine.mass, null).reduce_add();
                }
                else {
                    p.mass += select(mask, buffer.mass, null).reduce_add();
                }
            }
        }
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
