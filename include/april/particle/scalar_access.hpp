#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/source.hpp"

namespace april::particle::internal {
	template<ParticleField M, ParticleField N, IsParticleAttributes UserDataT> struct ScalarParticleRef;
	template <ParticleField ReadMask, ParticleField WriteMask> struct PackedParticleBuffer;

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
	// Ghost struct with references to particle data to abstract memory layout. Dissolves at compile time.
	template<ParticleField ReadMask, ParticleField WriteMask, IsParticleAttributes UserDataT>
    struct ScalarParticleRef {
    private:
       // Helper alias to pass the mutable type, the const type, and the field type
       template <typename MutT, typename ConstT, ParticleField F>
       using field_t = field_access_t<MutT, ConstT, F, ReadMask, WriteMask>;

       using MutVec3Ref   = math::Vec3Proxy<vec3::type>;
       using ConstVec3Ref = const math::Vec3Proxy<const vec3::type>;

    public:
		static constexpr ParticleField ReadAccess  = ReadMask;
		static constexpr ParticleField WriteAccess = WriteMask;

       // construct from ParticleSource
       template<class S>
       explicit ScalarParticleRef(const S & source) noexcept
          : force       (init_scalar_field<ReadMask, WriteMask, ParticleField::force>        (source))
          , position    (init_scalar_field<ReadMask, WriteMask, ParticleField::position>      (source))
          , velocity    (init_scalar_field<ReadMask, WriteMask, ParticleField::velocity>      (source))
          , old_position(init_scalar_field<ReadMask, WriteMask, ParticleField::old_position>   (source))
          , mass        (init_scalar_field<ReadMask, WriteMask, ParticleField::mass>         (source))
          , state       (init_scalar_field<ReadMask, WriteMask, ParticleField::state>        (source))
          , type        (init_scalar_field<ReadMask, WriteMask, ParticleField::type>         (source))
          , id          (init_scalar_field<ReadMask, WriteMask, ParticleField::id>           (source))
          , attributes  (init_scalar_field<ReadMask, WriteMask, ParticleField::attributes>   (source))
       {}

       // construct from a more permissive reference (e.g., converting a mutable ref to a view)
       template<ParticleField OtherReadMask, ParticleField OtherWriteMask>
       requires ((WriteMask & OtherWriteMask) == WriteMask) // Can only narrow write permissions, not expand
		&& ((OtherWriteMask | OtherReadMask) == (ReadMask | WriteMask))// must have the exact same fields
       explicit ScalarParticleRef(const ScalarParticleRef<OtherReadMask, OtherWriteMask, UserDataT>& r) noexcept
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

       // convenience method to drop all write permissions
       auto to_view() const noexcept {
          return ScalarParticleRef<ReadMask | WriteMask, ParticleField::none, UserDataT>(*this);
       }

		auto broadcast() const noexcept {
	       return PackedParticleBuffer<ReadMask, WriteMask>(*this);
       }

       // Data Fields: Mutable Type, Const Type, Field Enum
       AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::force>        force;
       AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::position>     position;
       AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::velocity>     velocity;
       AP_NO_UNIQUE_ADDRESS field_t<MutVec3Ref,     ConstVec3Ref,         ParticleField::old_position> old_position;

       AP_NO_UNIQUE_ADDRESS field_t<double&,        const double&,        ParticleField::mass>         mass;
       AP_NO_UNIQUE_ADDRESS field_t<ParticleState&, const ParticleState&, ParticleField::state>        state;
       AP_NO_UNIQUE_ADDRESS field_t<ParticleType&,  const ParticleType&,  ParticleField::type>         type;
       AP_NO_UNIQUE_ADDRESS field_t<ParticleID&,    const ParticleID&,    ParticleField::id>           id;
       AP_NO_UNIQUE_ADDRESS field_t<UserDataT&,     const UserDataT&,     ParticleField::attributes>   attributes;
    };




	//---------
	// CONCEPTS
	//---------
	template<typename T> struct is_scalar_particle_ref_impl : std::false_type {};

	template<ParticleField RM, ParticleField WM, typename U>
	struct is_scalar_particle_ref_impl<ScalarParticleRef<RM, WM, U>> : std::true_type {};
}


namespace april::particle {
	template<typename T>
	concept IsScalarParticleAccessor = internal::is_scalar_particle_ref_impl<std::remove_cvref_t<T>>::value;
}












