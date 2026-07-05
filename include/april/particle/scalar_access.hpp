/**
* @file scalar_access.hpp
 * @brief Scalar particle access abstraction via lightweight "ghost" references.
 *
 * ScalarParticleRef provides a unified AoS-like interface (p.position, p.velocity, etc.)
 * for kernels, regardless of whether the underlying storage is AoS, SoA, or AoSoA.
 *
 * It enforces Read/Write masks at compile time using the poisoning system from source.hpp.
 * Any access to an undeclared field results in a clear compile-time error.
 *
 * Because the struct consists entirely of references and proxies, it is completely
 * optimized away by the compiler, leaving only direct memory accesses.
 */

#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/source.hpp"
#include "april/particle/attributes.hpp"

namespace april::particle::internal {
	// forward declarations
	template<ParticleField M, ParticleField N, IsParticleAttributes UserDataT> struct ScalarParticleRef;
	template <ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes> struct PackedParticleBuffer;

	/**
	 * Helper to initialize field references from a data source.
	 * Maps valid fields to actual memory references and forbidden fields
	 * to the AccessForbidden poison type.
	 */
	template<ParticleField ReadMask, ParticleField WriteMask, ParticleField F, typename Source>
	constexpr decltype(auto) init_scalar_field(const Source& src) {
		if constexpr (particle::internal::has_field_v<ReadMask | WriteMask, F>) {
			return *src.template get<F>();
		} else {
			return internal::AccessForbidden<F>();
		}
	}


	//-------------------
	// PARTICLE REFERENCE
	//-------------------
	/**
	 * Lightweight proxy ("ghost") struct that gives kernels AoS-like scalar interface,
	 * regardless of actual memory layout.
	 *
	 * Member types are resolved at compile time:
	 *   - Mutable reference   if the field is in WriteMask
	 *   - Const reference     if the field is only in ReadMask
	 *   - AccessForbidden     if the field was not requested
	 *
	 * The id field is deliberately excluded from writing to prevent accidental mutation.
	 */
	template<ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes Attributes>
    struct ScalarParticleRef {
		static constexpr ParticleField ReadAccess  = ReadMask;

		// soft restriction so user does not need to zero out the id specifically when using ParticleField::all
		static constexpr ParticleField WriteAccess = WriteMask & ~ParticleField::id; // id is read-only
    private:
		//Maps the field to a mutable ref, const ref, or poison based on masks.
		template <typename MutT, typename ConstT, ParticleField F>
		using field_t = field_access_t<MutT, ConstT, F, ReadAccess, WriteAccess>;

		using MutVec3Ref   = math::Vec3Proxy<vec3::type>;
		using ConstVec3Ref = const math::Vec3Proxy<const vec3::type>;

    public:

		/**
		* Standard constructor from a ParticleSource (collection of raw pointers).
		*/
		template<class S>
		explicit ScalarParticleRef(const S & source) noexcept
		  : force       (init_scalar_field<ReadAccess, WriteAccess, ParticleField::force>        (source))
		  , position    (init_scalar_field<ReadAccess, WriteAccess, ParticleField::position>     (source))
		  , velocity    (init_scalar_field<ReadAccess, WriteAccess, ParticleField::velocity>     (source))
		  , old_position(init_scalar_field<ReadAccess, WriteAccess, ParticleField::old_position> (source))
		  , mass        (init_scalar_field<ReadAccess, WriteAccess, ParticleField::mass>         (source))
		  , state       (init_scalar_field<ReadAccess, WriteAccess, ParticleField::state>        (source))
		  , type        (init_scalar_field<ReadAccess, WriteAccess, ParticleField::type>         (source))
		  , id          (init_scalar_field<ReadAccess, WriteAccess, ParticleField::id>           (source))
		  , attributes  (init_scalar_field<ReadAccess, WriteAccess, ParticleField::attributes>   (source))
		{}

		/**
		* Narrowing constructor.
		* Allows converting a mutable reference to a read-only view.
		* Expansion of permissions is forbidden by the 'requires' clause.
		*/
		template<ParticleField OtherReadMask, ParticleField OtherWriteMask>
		requires ((WriteAccess & OtherWriteMask) == WriteAccess) // Can only narrow write permissions, not expand
		&& ((OtherWriteMask | OtherReadMask) == (ReadAccess | WriteAccess))// must have the exact same fields
		explicit ScalarParticleRef(const ScalarParticleRef<OtherReadMask, OtherWriteMask, Attributes>& r) noexcept
          : force       (r.force)
          , position    (r.position)
          , velocity    (r.velocity)
          , old_position(r.old_position)
          , mass        (r.mass)
          , state       (r.state)
          , type        (r.type)
          , id          (r.id)
          , attributes  (r.attributes)
		{}

		/**
		* Creates a read-only view of the current particle.
		*/
		auto to_view() const noexcept {
		  return ScalarParticleRef<ReadAccess | WriteAccess, ParticleField::none, Attributes>(*this);
		}

		/**
		* Broadcasts this scalar particle's data into a SIMD buffer.
		* Used in N x 1 interactions (Block vs Scalar).
		*/
		auto broadcast() const noexcept {
		   return PackedParticleBuffer<ReadAccess, WriteAccess, Attributes>(*this);
		}

		// Data Fields: Mutable Type, Const Type, Field Enum
		// APRIL_NO_UNIQUE_ADDRESS ensures forbidden fields do not increase object size.
		APRIL_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::force>        force;
		APRIL_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::position>     position;
		APRIL_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::velocity>     velocity;
		APRIL_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::old_position> old_position;

		APRIL_NO_UNIQUE_ADDRESS field_t<double&,        const double&,        ParticleField::mass>         mass;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleState&, const ParticleState&, ParticleField::state>        state;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleType&,  const ParticleType&,  ParticleField::type>         type;
		APRIL_NO_UNIQUE_ADDRESS field_t<ParticleID&,    const ParticleID&,    ParticleField::id>           id;
		APRIL_NO_UNIQUE_ADDRESS field_t<Attributes&,    const Attributes&,    ParticleField::attributes>   attributes;
    };



	//---------
	// CONCEPTS
	//---------
	template<typename T> struct is_scalar_particle_ref_impl : std::false_type {};

	template<ParticleField RM, ParticleField WM, typename U>
	struct is_scalar_particle_ref_impl<ScalarParticleRef<RM, WM, U>> : std::true_type {};
}


namespace april::particle {
	/**
	 * @brief Concept to identify any ScalarParticleRef regardless of its masks.
	 */
	template<typename T>
	concept IsScalarParticleAccessor = internal::is_scalar_particle_ref_impl<std::remove_cvref_t<T>>::value;
}


