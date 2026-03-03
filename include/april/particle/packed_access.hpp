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

            if constexpr (particle::internal::has_field_v<ReadMask, ParticleField::position>) {
                buf.position.x = scalar.position.x;
                buf.position.y = scalar.position.y;
                buf.position.z = scalar.position.z;
            } else if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::position>) {
                buf.position = pvec3(0.0);
            }

            if constexpr (particle::internal::has_field_v<ReadMask, ParticleField::old_position>) {
                buf.old_position.x = scalar.old_position.x;
                buf.old_position.y = scalar.old_position.y;
                buf.old_position.z = scalar.old_position.z;
            } else if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::old_position>) {
                buf.old_position = pvec3(0.0);
            }

            if constexpr (particle::internal::has_field_v<ReadMask, ParticleField::velocity>) {
                buf.velocity.x = scalar.velocity.x;
                buf.velocity.y = scalar.velocity.y;
                buf.velocity.z = scalar.velocity.z;
            } else if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::velocity>) {
                buf.velocity = pvec3(0.0);
            }

            if constexpr (particle::internal::has_field_v<ReadMask, ParticleField::force>) {
                buf.force.x = scalar.force.x;
                buf.force.y = scalar.force.y;
                buf.force.z = scalar.force.z;
            } else if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::force>) {
                buf.force = pvec3(0.0);
            }

            if constexpr (particle::internal::has_field_v<ReadMask, ParticleField::mass>) {
                buf.mass = scalar.mass;
            } else if constexpr (particle::internal::has_field_v<WriteMask, ParticleField::mass>) {
                buf.mass = 0.0;
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



    //-------------------
    // PACKED ACCUMULATOR
    //-------------------
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedAccumulator {
        static constexpr ParticleField RWMask = ReadMask & WriteMask;
        static constexpr ParticleField WOMask = WriteMask & ~ReadMask;

        PackedParticleBuffer<ReadMask, WriteMask> buffer;

        // Pristine and Acc only track fields that are Read + Write
        AP_NO_UNIQUE_ADDRESS PackedParticleBuffer<RWMask, ParticleField::none> pristine;
        AP_NO_UNIQUE_ADDRESS PackedParticleBuffer<RWMask, ParticleField::none> acc;

        // Construct directly from a packed reference
        template<typename attr>
        AP_FORCE_INLINE explicit PackedAccumulator(const PackedParticleRef<ReadMask, WriteMask, attr>& packed_ref)
            : buffer(packed_ref.load_buffer()) {
            reset();
        }

        // Construct from an existing buffer (necessary for duplicating in symmetric self-interaction)
        AP_FORCE_INLINE explicit PackedAccumulator(const PackedParticleBuffer<ReadMask, WriteMask>& buf)
            : buffer(buf) {
            reset();
        }

        template <typename ScalarAccessor>
           requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
       AP_FORCE_INLINE explicit PackedAccumulator(const ScalarAccessor& scalar)
            : buffer(PackedParticleBuffer<ReadMask, WriteMask>::broadcast(scalar)) {
            reset();
        }


        AP_FORCE_INLINE void reset() {
            if constexpr (has_field_v<RWMask, ParticleField::position>) acc.position = pvec3(0.0);
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) acc.old_position = pvec3(0.0);
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) acc.velocity = pvec3(0.0);
            if constexpr (has_field_v<RWMask, ParticleField::force>) acc.force = pvec3(0.0);
            if constexpr (has_field_v<RWMask, ParticleField::mass>) acc.mass = 0.0;
        }

        AP_FORCE_INLINE void save_pristine() {
            if constexpr (has_field_v<RWMask, ParticleField::position>) pristine.position = buffer.position;
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) pristine.old_position = buffer.old_position;
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) pristine.velocity = buffer.velocity;
            if constexpr (has_field_v<RWMask, ParticleField::force>) pristine.force = buffer.force;
            if constexpr (has_field_v<RWMask, ParticleField::mass>) pristine.mass = buffer.mass;
        }

        AP_FORCE_INLINE void extract_and_restore() {
            if constexpr (has_field_v<RWMask, ParticleField::position>) {
                acc.position += (buffer.position - pristine.position);
                buffer.position = pristine.position;
            }
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                acc.old_position += (buffer.old_position - pristine.old_position);
                buffer.old_position = pristine.old_position;
            }
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                acc.velocity += (buffer.velocity - pristine.velocity);
                buffer.velocity = pristine.velocity;
            }
            if constexpr (has_field_v<RWMask, ParticleField::force>) {
                acc.force += (buffer.force - pristine.force);
                buffer.force = pristine.force;
            }
            if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                acc.mass += (buffer.mass - pristine.mass);
                buffer.mass = pristine.mass;
            }
        }

        template<unsigned K = 1>
        AP_FORCE_INLINE void rotate_right() {
            buffer.template rotate_right<K>();
            acc.template rotate_right<K>();
        }

        template<unsigned K = 1>
        AP_FORCE_INLINE void rotate_left() {
            buffer.template rotate_left<K>();
            acc.template rotate_left<K>();
        }

        AP_FORCE_INLINE auto to_view() {
            save_pristine();
            return buffer.to_view();
        }

        // Unmasked Write-Back
        template<typename Attr>
        AP_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr> & packed_ref) {
            // Apply accumulators
            if constexpr (has_field_v<RWMask, ParticleField::position>) buffer.position += acc.position;
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) buffer.old_position += acc.old_position;
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) buffer.velocity += acc.velocity;
            if constexpr (has_field_v<RWMask, ParticleField::force>) buffer.force += acc.force;
            if constexpr (has_field_v<RWMask, ParticleField::mass>) buffer.mass += acc.mass;

            // Flush to memory
            packed_ref.update(buffer);
        }

        // Masked Write-Back (For SIMD partial tails)
        template<typename Attr>
        AP_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr> & packed_ref, const packed_mask & mask) {
            const packed null = 0.0;

            // Apply accumulators strictly to valid lanes
            if constexpr (has_field_v<RWMask, ParticleField::position>) {
                buffer.position.x += select(mask, acc.position.x, null);
                buffer.position.y += select(mask, acc.position.y, null);
                buffer.position.z += select(mask, acc.position.z, null);
            }
            if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                buffer.old_position.x += select(mask, acc.old_position.x, null);
                buffer.old_position.y += select(mask, acc.old_position.y, null);
                buffer.old_position.z += select(mask, acc.old_position.z, null);
            }
            if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                buffer.velocity.x += select(mask, acc.velocity.x, null);
                buffer.velocity.y += select(mask, acc.velocity.y, null);
                buffer.velocity.z += select(mask, acc.velocity.z, null);
            }
            if constexpr (has_field_v<RWMask, ParticleField::force>) {
                buffer.force.x += select(mask, acc.force.x, null);
                buffer.force.y += select(mask, acc.force.y, null);
                buffer.force.z += select(mask, acc.force.z, null);
            }
            if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                buffer.mass += select(mask, acc.mass, null);
            }

            // Mask Write-Only fields directly in the buffer so junk isn't written to dead lanes
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                buffer.position.x = select(mask, buffer.position.x, null);
                buffer.position.y = select(mask, buffer.position.y, null);
                buffer.position.z = select(mask, buffer.position.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                buffer.old_position.x = select(mask, buffer.old_position.x, null);
                buffer.old_position.y = select(mask, buffer.old_position.y, null);
                buffer.old_position.z = select(mask, buffer.old_position.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                buffer.velocity.x = select(mask, buffer.velocity.x, null);
                buffer.velocity.y = select(mask, buffer.velocity.y, null);
                buffer.velocity.z = select(mask, buffer.velocity.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                buffer.force.x = select(mask, buffer.force.x, null);
                buffer.force.y = select(mask, buffer.force.y, null);
                buffer.force.z = select(mask, buffer.force.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                buffer.mass = select(mask, buffer.mass, null);
            }

            // Flush to memory
            packed_ref.update(buffer);
        }

        // Unmasked Reduction
        template <typename ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p) const
        {
            if constexpr (has_field_v<WriteMask, ParticleField::position>) {
                if constexpr (has_field_v<RWMask, ParticleField::position>) {
                    p.position.x += acc.position.x.reduce_add();
                    p.position.y += acc.position.y.reduce_add();
                    p.position.z += acc.position.z.reduce_add();
                } else {
                    p.position.x += buffer.position.x.reduce_add();
                    p.position.y += buffer.position.y.reduce_add();
                    p.position.z += buffer.position.z.reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::old_position>) {
                if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                    p.old_position.x += acc.old_position.x.reduce_add();
                    p.old_position.y += acc.old_position.y.reduce_add();
                    p.old_position.z += acc.old_position.z.reduce_add();
                } else {
                    p.old_position.x += buffer.old_position.x.reduce_add();
                    p.old_position.y += buffer.old_position.y.reduce_add();
                    p.old_position.z += buffer.old_position.z.reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::velocity>) {
                if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                    p.velocity.x += acc.velocity.x.reduce_add();
                    p.velocity.y += acc.velocity.y.reduce_add();
                    p.velocity.z += acc.velocity.z.reduce_add();
                } else {
                    p.velocity.x += buffer.velocity.x.reduce_add();
                    p.velocity.y += buffer.velocity.y.reduce_add();
                    p.velocity.z += buffer.velocity.z.reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::force>) {
                if constexpr (has_field_v<RWMask, ParticleField::force>) {
                    p.force.x += acc.force.x.reduce_add();
                    p.force.y += acc.force.y.reduce_add();
                    p.force.z += acc.force.z.reduce_add();
                } else {
                    p.force.x += buffer.force.x.reduce_add();
                    p.force.y += buffer.force.y.reduce_add();
                    p.force.z += buffer.force.z.reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::mass>) {
                if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                    p.mass += acc.mass.reduce_add();
                } else {
                    p.mass += buffer.mass.reduce_add();
                }
            }
        }

        // Masked Reduction (For partial tail logic)
       template <typename ScalarAccessor, typename MaskT>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p, const MaskT& mask) const
        {
            const packed null = 0.0;

            if constexpr (has_field_v<WriteMask, ParticleField::position>) {
                if constexpr (has_field_v<RWMask, ParticleField::position>) {
                    p.position.x += select(mask, acc.position.x, null).reduce_add();
                    p.position.y += select(mask, acc.position.y, null).reduce_add();
                    p.position.z += select(mask, acc.position.z, null).reduce_add();
                } else {
                    p.position.x += select(mask, buffer.position.x, null).reduce_add();
                    p.position.y += select(mask, buffer.position.y, null).reduce_add();
                    p.position.z += select(mask, buffer.position.z, null).reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::old_position>) {
                if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
                    p.old_position.x += select(mask, acc.old_position.x, null).reduce_add();
                    p.old_position.y += select(mask, acc.old_position.y, null).reduce_add();
                    p.old_position.z += select(mask, acc.old_position.z, null).reduce_add();
                } else {
                    p.old_position.x += select(mask, buffer.old_position.x, null).reduce_add();
                    p.old_position.y += select(mask, buffer.old_position.y, null).reduce_add();
                    p.old_position.z += select(mask, buffer.old_position.z, null).reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::velocity>) {
                if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
                    p.velocity.x += select(mask, acc.velocity.x, null).reduce_add();
                    p.velocity.y += select(mask, acc.velocity.y, null).reduce_add();
                    p.velocity.z += select(mask, acc.velocity.z, null).reduce_add();
                } else {
                    p.velocity.x += select(mask, buffer.velocity.x, null).reduce_add();
                    p.velocity.y += select(mask, buffer.velocity.y, null).reduce_add();
                    p.velocity.z += select(mask, buffer.velocity.z, null).reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::force>) {
                if constexpr (has_field_v<RWMask, ParticleField::force>) {
                    p.force.x += select(mask, acc.force.x, null).reduce_add();
                    p.force.y += select(mask, acc.force.y, null).reduce_add();
                    p.force.z += select(mask, acc.force.z, null).reduce_add();
                } else {
                    p.force.x += select(mask, buffer.force.x, null).reduce_add();
                    p.force.y += select(mask, buffer.force.y, null).reduce_add();
                    p.force.z += select(mask, buffer.force.z, null).reduce_add();
                }
            }

            if constexpr (has_field_v<WriteMask, ParticleField::mass>) {
                if constexpr (has_field_v<RWMask, ParticleField::mass>) {
                    p.mass += select(mask, acc.mass, null).reduce_add();
                } else {
                    p.mass += select(mask, buffer.mass, null).reduce_add();
                }
            }
        }

        // Shift reciprocal forces before merging (Symmetric self-interaction only)
        // template<unsigned K>
        // AP_FORCE_INLINE void merge_shifted_left(const DeltaAccumulator& other) {
        //     if constexpr (has_field_v<RWMask, ParticleField::position>) {
        //         acc.position.x += other.acc.position.x.template rotate_left<K>();
        //         acc.position.y += other.acc.position.y.template rotate_left<K>();
        //         acc.position.z += other.acc.position.z.template rotate_left<K>();
        //     }
        //     if constexpr (has_field_v<RWMask, ParticleField::old_position>) {
        //         acc.old_position.x += other.acc.old_position.x.template rotate_left<K>();
        //         acc.old_position.y += other.acc.old_position.y.template rotate_left<K>();
        //         acc.old_position.z += other.acc.old_position.z.template rotate_left<K>();
        //     }
        //     if constexpr (has_field_v<RWMask, ParticleField::velocity>) {
        //         acc.velocity.x += other.acc.velocity.x.template rotate_left<K>();
        //         acc.velocity.y += other.acc.velocity.y.template rotate_left<K>();
        //         acc.velocity.z += other.acc.velocity.z.template rotate_left<K>();
        //     }
        //     if constexpr (has_field_v<RWMask, ParticleField::force>) {
        //         acc.force.x += other.acc.force.x.template rotate_left<K>();
        //         acc.force.y += other.acc.force.y.template rotate_left<K>();
        //         acc.force.z += other.acc.force.z.template rotate_left<K>();
        //     }
        // }
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
