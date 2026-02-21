#pragma once

#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/particle/particle_types.hpp"
#include "april/particle/fields.hpp"
#include "april/particle/source.hpp"

namespace april::env {
	template<ParticleField M, IsUserData UserDataT> struct ScalarParticleView;
	template<ParticleField M, IsUserData UserDataT> struct ScalarParticleRef;


	template<ParticleField M, ParticleField F, typename Source>
	constexpr decltype(auto) init_scalar_field(const Source& src) {
		if constexpr (has_field_v<M, F>) {
			return *src.template get<F>();
		} else {
			return internal::AccessForbidden<F>();
		}
	}



	//-------------------
	// PARTICLE REFERENCE
	//-------------------
	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<ParticleField M, IsUserData UserDataT>
	struct ScalarParticleRef {
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
			, user_data   (init_scalar_field<M, ParticleField::user_data>		(source))
		{}

		ScalarParticleView<M, UserDataT> to_view() noexcept {
			return ScalarParticleView<M, UserDataT>(*this);
		}

		using vec3ref = math::Vec3Proxy<vec3::type>;

		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::force, M> force;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::position, M> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::velocity, M> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<vec3ref, ParticleField::old_position, M> old_position;
		AP_NO_UNIQUE_ADDRESS field_type_t<double&, ParticleField::mass, M> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<ParticleState&, ParticleField::state, M> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type, M> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id, M> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<UserDataT&, ParticleField::user_data, M> user_data;
	};


	//------------------------
	// RESTRICTED PARTICLE REF
	//------------------------
	// Restricted reference allowing only the force field to be modified, used for fields.
	template<ParticleField M, IsUserData UserDataT>
	struct ScalarRestrictedParticleRef {

		explicit ScalarRestrictedParticleRef(const auto& source) noexcept
			: force       (init_scalar_field<M, ParticleField::force>			(source))
			, position    (init_scalar_field<M, ParticleField::position>		(source))
			, velocity    (init_scalar_field<M, ParticleField::velocity>		(source))
			, old_position(init_scalar_field<M, ParticleField::old_position>	(source))
			, mass        (init_scalar_field<M, ParticleField::mass>			(source))
			, state       (init_scalar_field<M, ParticleField::state>			(source))
			, type        (init_scalar_field<M, ParticleField::type>			(source))
			, id          (init_scalar_field<M, ParticleField::id>				(source))
			, user_data   (init_scalar_field<M, ParticleField::user_data>		(source))
		{}

		ScalarParticleView<M, UserDataT> to_view() noexcept {
			return ScalarParticleView<M, UserDataT>(*this);
		}

		using Vec3Ref = math::Vec3Proxy<vec3::type>;
		using ConstVec3Ref = math::Vec3Proxy<const vec3::type>;


		// everything by const reference except for force
		Vec3Ref force;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::position, M> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::velocity, M> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::old_position, M> old_position;
		AP_NO_UNIQUE_ADDRESS field_type_t<const double&, ParticleField::mass, M> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleState, ParticleField::state, M> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type, M> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id, M> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<const UserDataT&, ParticleField::user_data, M> user_data;
	};



	//--------------
	// PARTICLE VIEW
	//--------------
	// Immutable reference to particle data, intended for read-only access (e.g., monitors).
	template<ParticleField M, IsUserData UserDataT>
	struct ScalarParticleView {

		explicit ScalarParticleView(const auto & source) noexcept
			: force       (init_scalar_field<M, ParticleField::force>			(source))
			, position    (init_scalar_field<M, ParticleField::position>		(source))
			, velocity    (init_scalar_field<M, ParticleField::velocity>		(source))
			, old_position(init_scalar_field<M, ParticleField::old_position>	(source))
			, mass        (init_scalar_field<M, ParticleField::mass>			(source))
			, state       (init_scalar_field<M, ParticleField::state>			(source))
			, type        (init_scalar_field<M, ParticleField::type>			(source))
			, id          (init_scalar_field<M, ParticleField::id>				(source))
			, user_data   (init_scalar_field<M, ParticleField::user_data>		(source))
		{}

		template<typename RefT>
		requires (
			std::is_same_v<std::remove_cvref_t<RefT>, ScalarParticleRef<M, UserDataT>> ||
			std::is_same_v<std::remove_cvref_t<RefT>, ScalarRestrictedParticleRef<M, UserDataT>>
		)
		explicit ScalarParticleView(const RefT& r) noexcept
			: force        (r.force)     // Implicitly converts vec3& -> const vec3&
			, position     (r.position)  // Exact match (const vec3& -> const vec3&)
			, velocity     (r.velocity)
			, old_position (r.old_position)
			, mass         (r.mass)      // Exact match (const double -> const double)
			, state        (r.state)
			, type         (r.type)
			, id           (r.id)
			, user_data    (r.user_data)
			{}

		using ConstVec3Ref = math::Vec3Proxy<const vec3::type>;

		// everything by const reference
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::force, M> force;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::position, M> position;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::velocity, M> velocity;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ConstVec3Ref, ParticleField::old_position, M> old_position;
		AP_NO_UNIQUE_ADDRESS field_type_t<const double&, ParticleField::mass, M> mass;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleState, ParticleField::state, M> state;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleType, ParticleField::type, M> type;
		AP_NO_UNIQUE_ADDRESS field_type_t<const ParticleID, ParticleField::id, M> id;
		AP_NO_UNIQUE_ADDRESS field_type_t<const UserDataT&, ParticleField::user_data, M> user_data;
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

	// Concepts
	template<typename T>
	concept IsScalarRestrictedRef = is_restricted_ref_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleRef = is_particle_ref_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleView = is_particle_view_impl<std::remove_cvref_t<T>>::value;

	template<typename T>
	concept IsScalarParticleAccessor =
		IsScalarRestrictedRef<T> ||
		IsScalarParticleRef<T> ||
		IsScalarParticleView<T>;

}



