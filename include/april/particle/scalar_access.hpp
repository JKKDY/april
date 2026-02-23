#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/source.hpp"

namespace april::particle::internal {
	template<ParticleField M, IsParticleAttributes UserDataT> struct ScalarParticleView;
	template<ParticleField M, IsParticleAttributes UserDataT> struct ScalarParticleRef;


	template<ParticleField M, ParticleField F, typename Source>
	constexpr decltype(auto) init_scalar_field(const Source& src) {
		if constexpr (particle::internal::has_field_v<M, F>) {
			return *src.template get<F>();
		} else {
			return internal::AccessForbidden<F>();
		}
	}



	//-------------------
	// PARTICLE REFERENCE
	//-------------------
	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<ParticleField M, IsParticleAttributes UserDataT>
	struct ScalarParticleRef {
	private:
		template <typename T, ParticleField F> using field_type_t = field_type_t<T, F, M>;
		using vec3ref = math::Vec3Proxy<vec3::type>;
	public:

		template<class S>
		explicit ScalarParticleRef(const S & source) noexcept
			: force       (init_scalar_field<M, ParticleField::force>			(source))
			, position    (init_scalar_field<M, ParticleField::position>		(source))
			, velocity    (init_scalar_field<M, ParticleField::velocity>		(source))
			, old_position(init_scalar_field<M, ParticleField::old_position>	(source))
			, mass        (init_scalar_field<M, ParticleField::mass>			(source))
			, state       (init_scalar_field<M, ParticleField::state>			(source))
			, type        (init_scalar_field<M, ParticleField::type>			(source))
			, id          (init_scalar_field<M, ParticleField::id>				(source))
			, attributes   (init_scalar_field<M, ParticleField::attributes>		(source))
		{}

		ScalarParticleView<M, UserDataT> to_view() noexcept {
			return ScalarParticleView<M, UserDataT>(*this);
		}

		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::force> force;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::position> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::velocity> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::old_position> old_position;
		AP_NO_UNIQUE_ADDRESS field_type_t<double&, ParticleField::mass> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<ParticleState&, ParticleField::state> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<UserDataT&, ParticleField::attributes> attributes;
	};


	//------------------------
	// RESTRICTED PARTICLE REF
	//------------------------
	// Restricted reference allowing only the force field to be modified, used for fields.
	template<ParticleField M, IsParticleAttributes UserDataT>
	struct ScalarRestrictedParticleRef {
	private:
	    template <typename T, ParticleField F> using field_type_t = field_type_t<T, F, M>;
	    using Vec3Ref = math::Vec3Proxy<vec3::type>;
	    using ConstVec3Ref = math::Vec3Proxy<const vec3::type>;
	public:

	    explicit ScalarRestrictedParticleRef(const auto& source) noexcept
	       : force       (init_scalar_field<M, ParticleField::force>        (source))
	       , position    (init_scalar_field<M, ParticleField::position>      (source))
	       , velocity    (init_scalar_field<M, ParticleField::velocity>      (source))
	       , old_position(init_scalar_field<M, ParticleField::old_position>   (source))
	       , mass        (init_scalar_field<M, ParticleField::mass>         (source))
	       , state       (init_scalar_field<M, ParticleField::state>        (source))
	       , type        (init_scalar_field<M, ParticleField::type>         (source))
	       , id          (init_scalar_field<M, ParticleField::id>          (source))
	       , attributes   (init_scalar_field<M, ParticleField::attributes>     (source))
	    {}

	    ScalarParticleView<M, UserDataT> to_view() noexcept {
	       return ScalarParticleView<M, UserDataT>(*this);
	    }

	    AP_NO_UNIQUE_ADDRESS field_type_t<Vec3Ref, ParticleField::force> force;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::position> position;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::velocity> velocity;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::old_position> old_position;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const double&, ParticleField::mass> mass;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleState, ParticleField::state> state;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type> type;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id> id;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const UserDataT&, ParticleField::attributes> attributes;
	};



	//--------------
	// PARTICLE VIEW
	//--------------
	// Immutable reference to particle data, intended for read-only access (e.g., monitors).
	template<ParticleField M, IsParticleAttributes UserDataT>
	struct ScalarParticleView {
	private:
	    template <typename T, ParticleField F> using field_type_t = field_type_t<T, F, M>;
	    using ConstVec3Ref = math::Vec3Proxy<const vec3::type>;
	public:

	    explicit ScalarParticleView(const auto & source) noexcept
	       : force       (init_scalar_field<M, ParticleField::force>        (source))
	       , position    (init_scalar_field<M, ParticleField::position>      (source))
	       , velocity    (init_scalar_field<M, ParticleField::velocity>      (source))
	       , old_position(init_scalar_field<M, ParticleField::old_position>   (source))
	       , mass        (init_scalar_field<M, ParticleField::mass>         (source))
	       , state       (init_scalar_field<M, ParticleField::state>        (source))
	       , type        (init_scalar_field<M, ParticleField::type>         (source))
	       , id          (init_scalar_field<M, ParticleField::id>          (source))
	       , attributes   (init_scalar_field<M, ParticleField::attributes>     (source))
	    {}

	    template<typename RefT>
	    requires (
	       std::is_same_v<std::remove_cvref_t<RefT>, ScalarParticleRef<M, UserDataT>> ||
	       std::is_same_v<std::remove_cvref_t<RefT>, ScalarRestrictedParticleRef<M, UserDataT>>
	    )
	    explicit ScalarParticleView(const RefT& r) noexcept
	       : force        (r.force)
	       , position     (r.position)
	       , velocity     (r.velocity)
	       , old_position (r.old_position)
	       , mass         (r.mass)
	       , state        (r.state)
	       , type         (r.type)
	       , id           (r.id)
	       , attributes    (r.attributes)
	       {}

	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::force> force;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::position> position;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::velocity> velocity;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::old_position> old_position;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const double&, ParticleField::mass> mass;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleState, ParticleField::state> state;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type> type;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id> id;
	    AP_NO_UNIQUE_ADDRESS field_type_t<const UserDataT&, ParticleField::attributes> attributes;
	};



	//---------
	// CONCEPTS
	//---------
	template<typename T> struct is_restricted_ref_impl : std::false_type {};
	template<typename T> struct is_particle_ref_impl   : std::false_type {};
	template<typename T> struct is_particle_view_impl  : std::false_type {};

	// Specialization for Accessors
	template<auto M, typename U> struct is_restricted_ref_impl<ScalarRestrictedParticleRef<M, U>> : std::true_type {};
	template<auto M, typename U> struct is_particle_ref_impl<ScalarParticleRef<M, U>> : std::true_type {};
	template<auto M, typename U> struct is_particle_view_impl<ScalarParticleView<M, U>> : std::true_type {};
}


namespace april::particle {
	// Concepts
	template<typename T>
	concept IsScalarRestrictedRef = internal::is_restricted_ref_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleRef = internal::is_particle_ref_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleView = internal::is_particle_view_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleAccessor =
		IsScalarRestrictedRef<T> ||
		IsScalarParticleRef<T> ||
		IsScalarParticleView<T>;
}












