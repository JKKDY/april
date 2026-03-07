#pragma once
#include "april/base/types.hpp"
#include "april/particle/source.hpp"
#include "april/simd/packed.hpp"
#include "april/simd/packed_ref.hpp"
#include "april/math/vec3.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/scalar_access.hpp"
#include "april/particle/attributes.hpp"

namespace april::particle::internal {
    // fwd declaration
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedBufferView;
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes U>
    struct PackedParticleRef;


    //--------------------
    // PACKED PARTICLE REF
    //--------------------
    // holds packed references which act like packed types but will touch memory on loads/writebacks.
    // for better perf use buffers for a single load at the beginning and single write back at the end
    template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
    struct PackedParticleRef {
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask;
    private:
        // helper for switching between data pointers and poison
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

        //switch between mutable and const type depending on read and write mask
        template <typename MutT, typename ConstT, ParticleField F>
        using field_t = field_access_t<MutT, ConstT, F, ReadAccess, WriteAccess>;

        template <typename T>
        using target_reg_t = std::conditional_t<
            std::is_floating_point_v<T>,
            packed::value_type,
            std::conditional_t<std::is_signed_v<T>, packedi::value_type, packedu::value_type>
        >;

        // declare a packed type
        template <typename MemT, ParticleField F>
        using packed_field_t = field_t<
            simd::PackedRef<MemT, simd::Packed<target_reg_t<MemT>>>,
            const simd::PackedRef<const MemT, simd::Packed<target_reg_t<MemT>>>,
            F
        >;

        // declare a vec type
        template <ParticleField F>
        using vec3_field_t = field_t<math::Vec3Proxy<pvec3::type>, const math::Vec3Proxy<const pvec3::type>, F>;

        // declare a raw pointer
        template<typename T, ParticleField F>
        using Ptr = field_access_t<T* AP_RESTRICT, const T* AP_RESTRICT, F, ReadAccess, WriteAccess>;

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
            {}

        // Cross-constructor for narrowing write permissions (e.g. creating a view)
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
        {}

        // a view is just a PackedParticleRef with no write permissions
        auto to_view() const noexcept {
            return PackedParticleRef<ReadAccess | WriteAccess, ParticleField::none, Attributes>(*this);
        }

        PackedParticleBuffer<ReadAccess, WriteAccess> load_buffer() const noexcept {
            return PackedParticleBuffer<ReadAccess, WriteAccess>(*this);
        }

        // Data members with strict const-correctness
        AP_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::force> force;
        AP_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS vec3_field_t<ParticleField::old_position> old_position;

        AP_NO_UNIQUE_ADDRESS packed_field_t<double,        ParticleField::mass> mass;
        AP_NO_UNIQUE_ADDRESS packed_field_t<std::underlying_type_t<ParticleState>, ParticleField::state> state;
        AP_NO_UNIQUE_ADDRESS packed_field_t<ParticleType,  ParticleField::type> type;
        AP_NO_UNIQUE_ADDRESS packed_field_t<ParticleID,    ParticleField::id> id;

        AP_NO_UNIQUE_ADDRESS Ptr<Attributes,    ParticleField::attributes> attributes;
    };


    //-----------------------
    // PACKED PARTICLE BUFFER
    //-----------------------
    // shadow object with actual SIMD registers. Allows for manipulations without direct write backs
    template <ParticleField ReadMask, ParticleField WriteMask>
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
    public:
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask;

        static constexpr ParticleField RWMask = ReadMask & WriteMask;  // read & write
        static constexpr ParticleField WOMask = WriteMask & ~ReadMask; // write only
        static constexpr ParticleField ROMask = ReadMask & ~WriteMask; // read only

        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::position> position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::old_position> old_position;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::velocity> velocity;
        AP_NO_UNIQUE_ADDRESS pvec3_t<ParticleField::force> force;

        AP_NO_UNIQUE_ADDRESS packed_float_t<ParticleField::mass> mass;
        AP_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::state> state;
        AP_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::type>  type;
        AP_NO_UNIQUE_ADDRESS packed_int_t<ParticleField::id>    id;

        PackedParticleBuffer() = default;

        // if in read mask read from memory else if (only) in write mask zero initialize
        template <typename attr>
        explicit PackedParticleBuffer(const PackedParticleRef<ReadMask, WriteMask, attr>& source) {
            if constexpr (has_field_v<ReadMask, ParticleField::position>) position = source.position;
            if constexpr (has_field_v<ReadMask, ParticleField::old_position>) old_position = source.old_position;
            if constexpr (has_field_v<ReadMask, ParticleField::velocity>) velocity = source.velocity;
            if constexpr (has_field_v<ReadMask, ParticleField::force>) force = source.force;
            if constexpr (has_field_v<ReadMask, ParticleField::mass>) mass = source.mass;
            if constexpr (has_field_v<ReadMask, ParticleField::state>) state = source.state;
            if constexpr (has_field_v<ReadMask, ParticleField::type>) type = source.type;
            if constexpr (has_field_v<ReadMask, ParticleField::state>) id = source.id;

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
            } else if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                old_position = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::velocity>) {
                velocity.x = scalar.velocity.x;
                velocity.y = scalar.velocity.y;
                velocity.z = scalar.velocity.z;
            } else if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                velocity = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::force>) {
                force.x = scalar.force.x;
                force.y = scalar.force.y;
                force.z = scalar.force.z;
            } else if constexpr (has_field_v<WOMask, ParticleField::force>) {
                force = pvec3(0.0);
            }

            if constexpr (has_field_v<ReadMask, ParticleField::mass>) {
                mass = scalar.mass;
            } else if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass = 0.0;
            }

            if constexpr (has_field_v<ReadMask, ParticleField::state>) {
                state = scalar.state;
            }

            if constexpr (has_field_v<ReadMask, ParticleField::type>) {
                type = scalar.type;
            }

            if constexpr (has_field_v<ReadMask, ParticleField::id>) {
                id = scalar.id;
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
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::state>) {
                state = state.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::type>) {
                type = type.template rotate_left<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::id>) {
                id = id.template rotate_left<K>();
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
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::type>) {
                type = type.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::state>) {
                state = state.template rotate_right<K>();
            }
            if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, ParticleField::id>) {
                id = id.template rotate_right<K>();
            }
        }

        // Accumulate reciprocal deltas from another buffer (Write-Only fields strictly)
        AP_FORCE_INLINE void accumulate(const PackedParticleBuffer& other) {
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                position += other.position;
            }
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                old_position += other.old_position;
            }
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                velocity += other.velocity;
            }
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                force += other.force;
            }
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass += other.mass;
            }
        }

        // Masked accumulation
        template <typename MaskT>
        AP_FORCE_INLINE void accumulate(const PackedParticleBuffer& other, const MaskT& mask) {
            const packed null = 0.0;
            if constexpr (has_field_v<WOMask, ParticleField::position>) {
                position.x += select(mask, other.position.x, null);
                position.y += select(mask, other.position.y, null);
                position.z += select(mask, other.position.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::old_position>) {
                old_position.x += select(mask, other.old_position.x, null);
                old_position.y += select(mask, other.old_position.y, null);
                old_position.z += select(mask, other.old_position.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::velocity>) {
                velocity.x += select(mask, other.velocity.x, null);
                velocity.y += select(mask, other.velocity.y, null);
                velocity.z += select(mask, other.velocity.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::force>) {
                force.x += select(mask, other.force.x, null);
                force.y += select(mask, other.force.y, null);
                force.z += select(mask, other.force.z, null);
            }
            if constexpr (has_field_v<WOMask, ParticleField::mass>) {
                mass += select(mask, other.mass, null);
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

            // STATE
            if constexpr (has_field_v<RWMask, ParticleField::state>) {
                packed_ref.state = state;
            }

            // TYPE
            if constexpr (has_field_v<RWMask, ParticleField::state>) {
                packed_ref.type = type;
            }

            // id is not assignable
        }

      // Masked SIMD Write-Back (For Tail Chunks)
        template <typename Attr, typename MaskT>
        AP_FORCE_INLINE void update_into(PackedParticleRef<ReadMask, WriteMask, Attr>& packed_ref, const MaskT & mask) const {
            update_vec_masked<ParticleField::position>(packed_ref.position, position, mask);
            update_vec_masked<ParticleField::old_position>(packed_ref.old_position, old_position, mask);
            update_vec_masked<ParticleField::velocity>(packed_ref.velocity, velocity, mask);
            update_vec_masked<ParticleField::force>(packed_ref.force, force, mask);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                packed_ref.mass += select(mask, mass, 0.0);
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                packed_ref.mass = select(mask, mass, packed_ref.mass);
        }

        // Unmasked Scalar Reduction (For Broadcast Buffers)
        template <typename ScalarAccessor>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p) const {
            reduce_vec_unmasked<ParticleField::position>(p.position, position);
            reduce_vec_unmasked<ParticleField::old_position>(p.old_position, old_position);
            reduce_vec_unmasked<ParticleField::velocity>(p.velocity, velocity);
            reduce_vec_unmasked<ParticleField::force>(p.force, force);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                p.mass += mass.reduce_add();
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot reduce RW mass.");
        }

        // Masked Scalar Reduction (For partial tail logic vs broadcasted scalars)
        template <typename ScalarAccessor, typename MaskT>
        requires april::particle::IsScalarParticleAccessor<ScalarAccessor>
        AP_FORCE_INLINE void reduce_into(ScalarAccessor& p, const MaskT& mask) const {
            reduce_vec_masked<ParticleField::position>(p.position, position, mask);
            reduce_vec_masked<ParticleField::old_position>(p.old_position, old_position, mask);
            reduce_vec_masked<ParticleField::velocity>(p.velocity, velocity, mask);
            reduce_vec_masked<ParticleField::force>(p.force, force, mask);

            // MASS (Scalar inline)
            if constexpr (has_field_v<WOMask, ParticleField::mass>)
                p.mass += select(mask, mass, 0.0).reduce_add();
            else if constexpr (has_field_v<RWMask, ParticleField::mass>)
                static_assert(sizeof(ScalarAccessor) == 0, "FATAL: Cannot masked reduce RW mass.");
        }

    private:
        // Unified Vector Write-Back (Masked)
        template <ParticleField F, typename DestT, typename SrcT, typename MaskT>
        AP_FORCE_INLINE static void update_vec_masked(DestT& dest, const SrcT& src, const MaskT& mask) {
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

        // Unified Vector Reduce (Unmasked)
        template <ParticleField F, typename ScalarT, typename SimdT>
        AP_FORCE_INLINE static void reduce_vec_unmasked(ScalarT& dest, const SimdT& src) {
            if constexpr (has_field_v<WOMask, F>) {
                dest.x += src.x.reduce_add();
                dest.y += src.y.reduce_add();
                dest.z += src.z.reduce_add();
            } else if constexpr (has_field_v<RWMask, F>) {
                static_assert(sizeof(ScalarT) == 0, "FATAL: Cannot reduce a Read-Write vector field from a SIMD register to a scalar.");
            }
        }

        // Unified Vector Reduce (Masked)
        template <ParticleField F, typename ScalarT, typename SimdT, typename MaskT>
        AP_FORCE_INLINE static void reduce_vec_masked(ScalarT& dest, const SimdT& src, const MaskT& mask) {
            if constexpr (has_field_v<WOMask, F>) {
                const packed null = 0.0;
                dest.x += select(mask, src.x, null).reduce_add();
                dest.y += select(mask, src.y, null).reduce_add();
                dest.z += select(mask, src.z, null).reduce_add();
            } else if constexpr (has_field_v<RWMask, F>) {
                static_assert(sizeof(ScalarT) == 0, "FATAL: Cannot perform masked reduction on a Read-Write vector field.");
            }
        }
    };


    //------------
    // BUFFER VIEW
    //------------
    // Views enforce the read-write rules on buffers
    template <ParticleField ReadMask, ParticleField WriteMask>
    struct PackedBufferView {
    private:
        template <ParticleField F, typename T>
          using view_ref_t = std::conditional_t<
              std::is_same_v<T, AccessForbidden<F>>, // if it's poison, keep it as poison (by value)
              T,
              std::conditional_t<has_field_v<WriteMask, F>, T&, const T&>
              // if it's valid and in WriteMask, make it a mutable ref else a const ref
          >;

        template <ParticleField F, typename T>
        using ref_t = view_ref_t<F, T>;

        using Buffer = PackedParticleBuffer<ReadMask, WriteMask>;

    public:
        static constexpr ParticleField ReadAccess  = ReadMask;
        static constexpr ParticleField WriteAccess = WriteMask;

        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::position, decltype(Buffer::position)> position;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::old_position, decltype(Buffer::old_position)> old_position;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::velocity, decltype(Buffer::velocity)> velocity;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::force, decltype(Buffer::force)> force;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::mass, decltype(Buffer::mass)> mass;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::state, decltype(Buffer::state)> state;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::type, decltype(Buffer::type)> type;
        AP_NO_UNIQUE_ADDRESS ref_t<ParticleField::id, decltype(Buffer::id)> id;

        // Directly maps buffer fields to references or poison copies. Zero branching.
        AP_FORCE_INLINE explicit PackedBufferView(Buffer& buf)
            : position(buf.position),
              old_position(buf.old_position),
              velocity(buf.velocity),
              force(buf.force),
              mass(buf.mass),
              state(buf.state),
              type(buf.type),
              id(buf.id)
            {}
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

    static_assert(IsPackedParticleRef<internal::PackedParticleRef<ParticleField::all, ParticleField::all, NoParticleAttributes>>);
    static_assert(IsPackedParticleView<internal::PackedBufferView<ParticleField::all, ParticleField::all>>);
    static_assert(IsPackedParticleBuffer<internal::PackedParticleBuffer<ParticleField::all, ParticleField::all>>);
}
