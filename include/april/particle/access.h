#pragma once

#include <string>
#include <sstream>

#include "april/common.h"
#include "april/particle/defs.h"
#include "april/particle/fields.h"

namespace april::env {

	//---------------------
	// PARTICLE DATA SOURCE
	//---------------------


	// field attribute helpers
	template<typename T, Field F, FieldMask M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, std::monostate>;

	template<FieldMask M, IsUserData U, bool IsConst>
	struct ParticleSource {

	    // selects T* or const T*
	    template<typename T>
	    using Ptr = std::conditional_t<IsConst, const T*, T*>;

	    // data pointers (optimized away if not in M)
	    [[no_unique_address]] field_type_t<Ptr<vec3>, Field::force, M> force;
	    [[no_unique_address]] field_type_t<Ptr<vec3>, Field::position, M> position;
	    [[no_unique_address]] field_type_t<Ptr<vec3>, Field::velocity, M> velocity;
	    [[no_unique_address]] field_type_t<Ptr<vec3>, Field::old_position, M> old_position;
	    [[no_unique_address]] field_type_t<Ptr<double>, Field::mass, M> mass;
	    [[no_unique_address]] field_type_t<Ptr<ParticleState>, Field::state, M> state;
	    [[no_unique_address]] field_type_t<Ptr<ParticleType>, Field::type, M> type;
	    [[no_unique_address]] field_type_t<Ptr<ParticleID>, Field::id, M> id;
	    [[no_unique_address]] field_type_t<Ptr<U>, Field::user_data, M> user_data;

		// getter (Used by ParticleRef/View)
		template<Field F>
		constexpr auto get() const {
			if constexpr (has_field_v<M, F>) {
				if constexpr (F == Field::force) return force;
				else if constexpr (F == Field::position) return position;
				else if constexpr (F == Field::velocity) return velocity;
				else if constexpr (F == Field::old_position) return old_position;
				else if constexpr (F == Field::mass) return mass;
				else if constexpr (F == Field::state) return state;
				else if constexpr (F == Field::type) return type;
				else if constexpr (F == Field::id) return id;
				else if constexpr (F == Field::user_data) return user_data;
			} else {
				return static_cast<Ptr<void>>(nullptr);
			}
		}
	};

	template<FieldMask M, Field F, typename Source>
	constexpr decltype(auto) init_field(const Source& src) {
		if constexpr (has_field_v<M, F>) {
			return *src.template get<F>();
		} else {
			return std::monostate{};
		}
	}


	//---------------------
	// CONTROLLED ACCESSORS
	//---------------------

	// forward declaration
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView;

	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef;


	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef {
		template<class S>
		explicit ParticleRef(const S & source)
			: force       (init_field<M, Field::force>			(source))
			, position    (init_field<M, Field::position>		(source))
			, velocity    (init_field<M, Field::velocity>		(source))
			, old_position(init_field<M, Field::old_position>	(source))
			, mass        (init_field<M, Field::mass>			(source))
			, state       (init_field<M, Field::state>			(source))
			, type        (init_field<M, Field::type>			(source))
			, id          (init_field<M, Field::id>				(source))
			, user_data   (init_field<M, Field::user_data>		(source))
		{}

		ParticleView<M, UserDataT> to_view() {
			return ParticleView<M, UserDataT>(*this);
		}

		[[no_unique_address]] field_type_t<vec3&, Field::force, M> force;
		[[no_unique_address]] field_type_t<vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<ParticleState&, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<UserDataT&, Field::user_data, M> user_data;
	};

	// Restricted reference allowing only the force field to be modified, used for fields.
	template<FieldMask M, IsUserData UserDataT>
	struct RestrictedParticleRef {

		explicit RestrictedParticleRef(const auto& source)
			: force       (init_field<M, Field::force>			(source))
			, position    (init_field<M, Field::position>		(source))
			, velocity    (init_field<M, Field::velocity>		(source))
			, old_position(init_field<M, Field::old_position>	(source))
			, mass        (init_field<M, Field::mass>			(source))
			, state       (init_field<M, Field::state>			(source))
			, type        (init_field<M, Field::type>			(source))
			, id          (init_field<M, Field::id>				(source))
			, user_data   (init_field<M, Field::user_data>		(source))
		{}

		ParticleView<M, UserDataT> to_view() {
			return ParticleView<M, UserDataT>(*this);
		}


		// everything by const reference except for force
		vec3& force;
		[[no_unique_address]] field_type_t<const vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<const vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<const double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<const ParticleState, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<const UserDataT&, Field::user_data, M> user_data;
	};


	// Immutable reference to particle data, intended for read-only access (e.g., monitors).
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView {

		explicit ParticleView(const auto & source)
			: force       (init_field<M, Field::force>			(source))
			, position    (init_field<M, Field::position>		(source))
			, velocity    (init_field<M, Field::velocity>		(source))
			, old_position(init_field<M, Field::old_position>	(source))
			, mass        (init_field<M, Field::mass>			(source))
			, state       (init_field<M, Field::state>			(source))
			, type        (init_field<M, Field::type>			(source))
			, id          (init_field<M, Field::id>				(source))
			, user_data   (init_field<M, Field::user_data>		(source))
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