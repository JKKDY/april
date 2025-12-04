#pragma once

#include <string>
#include <sstream>
#include <utility>

#include "april/common.h"
#include "april/particle/defs.h"
#include "april/particle/fields.h"

namespace april::env {

	//------------------
	// COMPILE-TIME DEFS
	//------------------

	// Fetcher concepts
	template<typename F>
	concept IsFetcher =
	requires { typename std::remove_reference_t<F>::UserDataT; } &&
	IsUserData<typename std::remove_reference_t<F>::UserDataT> &&
	requires(F f) {
		{ f.force() }         -> std::same_as<vec3&>;
		{ f.position() }      -> std::same_as<vec3&>;
		{ f.velocity() }      -> std::same_as<vec3&>;
		{ f.old_position() }  -> std::same_as<vec3&>;
		{ f.old_force() }     -> std::same_as<vec3&>;
		{ f.mass() }          -> std::same_as<double&>;
		{ f.state() }         -> std::same_as<ParticleState&>;
		{ f.type() }          -> std::same_as<ParticleType&>;
		{ f.id() }            -> std::same_as<ParticleID&>;
		{ f.user_data() }     -> std::same_as<typename std::remove_reference_t<F>::UserDataT&>;
	} &&
	std::is_copy_constructible_v<F> &&
	std::is_move_constructible_v<F>;


	// field attribute helpers
	template<typename T, Field F, FieldMask M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, std::monostate>;

	template< typename T, Field F, FieldMask M, typename Get>
	constexpr field_type_t<T, F, M> init_field(Get && get) noexcept(noexcept(get())) {
		if constexpr (has_field_v<M, F>) return std::forward<Get>(get)();  // call get
		else return std::monostate();
	}


	//---------------------
	// CONTROLLED ACCESSORS
	//---------------------

	// forward declaration
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView;

	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef {

		ParticleRef(const ParticleRef& other) = default;
		ParticleRef(ParticleRef&& other) = default;

		template<IsFetcher F>
		explicit ParticleRef(F && f)
			: force      ( init_field<vec3&,         		Field::force,		M>([&]() noexcept -> vec3&{ return f.force(); }) )
			, position   ( init_field<vec3&,         		Field::position,	M>([&]() noexcept -> vec3&{ return f.position(); }) )
			, velocity   ( init_field<vec3&,         		Field::velocity,	M>([&]() noexcept -> vec3&{ return f.velocity(); }) )
			, old_position(init_field<vec3&,         		Field::old_position,M>([&]() noexcept -> vec3&{ return f.old_position(); }) )
			, old_force  ( init_field<vec3&,         		Field::old_force, 	M>([&]() noexcept -> vec3&{ return f.old_force(); }) )
			, mass       ( init_field<double&,       		Field::mass,      	M>([&]() noexcept -> double&{ return f.mass(); }) )
			, state      ( init_field<ParticleState&,		Field::state,     	M>([&]() noexcept -> ParticleState&{ return f.state(); }) )
			, type       ( init_field<const ParticleType, 	Field::type,		M>([&]() noexcept -> ParticleType{ return f.type(); }) )
			, id         ( init_field<const ParticleID,   	Field::id,			M>([&]() noexcept -> ParticleID{ return f.id(); }) )
			, user_data  ( init_field<UserDataT&,			Field::user_data,	M>([&]() noexcept -> UserDataT&{ return f.user_data(); }) )
		{}

		ParticleView<M, UserDataT> to_view() {
			return ParticleView<M, UserDataT>(*this);
		}

		// everything by reference. mutable.
		[[no_unique_address]] field_type_t<vec3&, Field::force, M> force;
		[[no_unique_address]] field_type_t<vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<ParticleState&, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<UserDataT&, Field::user_data, M> user_data;
	};

	// Restricted reference allowing only the force field to be modified, used for fields.
	template<FieldMask M, IsUserData UserDataT>
	struct RestrictedParticleRef {

		template<IsFetcher F>
		explicit RestrictedParticleRef(F && f) // forwarding reference
			: force        ( f.force() )
			, position     ( init_field<const vec3&,		Field::position,	M>([&]() noexcept -> const vec3&{ return f.position(); }) )
			, velocity     ( init_field<const vec3&,		Field::velocity,	M>([&]() noexcept -> const vec3&{ return f.velocity(); }) )
			, old_position ( init_field<const vec3&,		Field::old_position,M>([&]() noexcept -> const vec3&{ return f.old_position(); }) )
			, old_force    ( init_field<const vec3&,		Field::old_force,	M>([&]() noexcept -> const vec3&{ return f.old_force(); }) )
			, mass         ( init_field<const double&,		Field::mass,		M>([&]() noexcept -> double { return f.mass(); }) )
			, state        ( init_field<const ParticleState,Field::state,		M>([&]() noexcept -> ParticleState { return f.state(); }) )
			, type         ( init_field<const ParticleType, Field::type,		M>([&]() noexcept -> ParticleType { return f.type(); }) )
			, id           ( init_field<const ParticleID,	Field::id,			M>([&]() noexcept -> ParticleID { return f.id(); }) )
			, user_data    ( init_field<const UserDataT&,	Field::user_data,	M>([&]() noexcept -> const UserDataT& { return f.user_data(); }) ) {
			static_assert(has_field_v<M, Field::force>,
			"ParticleRef must have Field::force to be converted to RestrictedParticleRef. "
			"Have you checked if all force fields have Field::force set?");
		}

		ParticleView<M, UserDataT> to_view() {
			return ParticleView<M, UserDataT>(*this);
		}

		// everything by const reference except for force
		vec3& force;
		[[no_unique_address]] field_type_t<const vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<const vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<const double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<const ParticleState, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<const UserDataT&, Field::user_data, M> user_data;
	};


		// Immutable reference to particle data, intended for read-only access (e.g., monitors).
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView {

		template<IsFetcher F>
		explicit ParticleView(F && f)
		   : force        ( init_field<const vec3&, 		Field::force,		M>([&]() noexcept -> const vec3&{ return f.force(); }) )
		   , position     ( init_field<const vec3&, 		Field::position,	M>([&]() noexcept -> const vec3&{ return f.position(); }) )
		   , velocity     ( init_field<const vec3&, 		Field::velocity,	M>([&]() noexcept -> const vec3&{ return f.velocity(); }) )
		   , old_position ( init_field<const vec3&, 		Field::old_position,M>([&]() noexcept -> const vec3&{ return f.old_position(); }) )
		   , old_force    ( init_field<const vec3&, 		Field::old_force,	M>([&]() noexcept -> const vec3&{ return f.old_force(); }) )
		   , mass         ( init_field<const double&,		Field::mass,		M>([&] noexcept { return f.mass(); }) )
		   , state        ( init_field<const ParticleState, Field::state,		M>([&] noexcept { return f.state(); }) )
		   , type         ( init_field<const ParticleType,	Field::type,		M>([&] noexcept { return f.type(); }) )
		   , id           ( init_field<const ParticleID,	Field::id,			M>([&] noexcept { return f.id(); }) )
		   , user_data    ( init_field<const UserDataT&,	Field::user_data,	M>([&]() noexcept -> const UserDataT&{ return f.user_data(); }) )
			{}

		template<typename RefT>
		requires (
			std::is_same_v<std::remove_cvref_t<RefT>, ParticleRef<M, UserDataT>> ||
			std::is_same_v<std::remove_cvref_t<RefT>, RestrictedParticleRef<M, UserDataT>>
		)
		explicit ParticleView(const RefT& r)
			: force        (r.force)     // Implicitly converts vec3& -> const vec3&
			, position     (r.position)  // Exact match (const vec3& -> const vec3&)
			, velocity     (r.velocity)
			, old_position (r.old_position)
			, old_force    (r.old_force)
			, mass         (r.mass)      // Exact match (const double -> const double)
			, state        (r.state)
			, type         (r.type)
			, id           (r.id)
			, user_data    (r.user_data)
				{}

		// everything by const reference
		[[no_unique_address]] field_type_t<const vec3&, Field::force, M> force;
		[[no_unique_address]] field_type_t<const vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<const vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<const double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<const ParticleState, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<const UserDataT&, Field::user_data, M> user_data;
	};


	// TODO maybe stick this somewhere else? like in particle.h
	// easy terminal diagnostics
	template<typename P>
	std::string particle_to_string(const P & p) {
		std::ostringstream oss;
		oss << "Particle ID: " << p.id << "\n"
			<< "Position: " << p.position.to_string() << "\n"
			<< "Velocity: " << p.velocity.to_string() << "\n"
			<< "Force: " << p.force.to_string() << "\n"
			<< "Mass: " << p.mass << "\n"
			<< "Type: " << p.type << "\n"
			<< "State: " << static_cast<int>(p.state) << "\n";
		return oss.str();
	}
}