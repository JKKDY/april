#pragma once

#include <string>
#include <type_traits>
#include <cstdint>
#include <sstream>
#include <utility>

#include "april/common.h"
#include "april/particle/defs.h"


namespace april::env {

	using FieldMask = uint32_t;

	enum class Field : FieldMask {
		none			= 0u,
		position     	= 1u << 0,
		velocity     	= 1u << 1,
		force        	= 1u << 2,
		old_position 	= 1u << 3,
		old_force    	= 1u << 4,
		state        	= 1u << 5,
		mass         	= 1u << 6,
		type         	= 1u << 7,
		id           	= 1u << 8,
		user_data    	= 1u << 9,
		all			 	= ~0u
	};

	constexpr FieldMask to_field_mask(Field f) { return static_cast<FieldMask>(f); }

	constexpr FieldMask operator|(const Field a, const Field b) { return to_field_mask(a) | to_field_mask(b); }
	constexpr FieldMask operator|(const FieldMask m, const Field f) { return m | to_field_mask(f); }
	constexpr FieldMask operator|(const Field a, const FieldMask m) { return to_field_mask(a) | m; }


	template<class T>
	concept HasFields = requires { std::remove_cvref_t<T>::fields; };

	template<HasFields Self>
	inline constexpr FieldMask FieldOf = std::remove_cvref_t<Self>::fields;

	template<FieldMask M, Field F>
	inline constexpr bool has_field_v = (M & to_field_mask(F)) != 0;




	template<typename T, Field F, FieldMask M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, std::monostate>;

	template< typename T, Field F, FieldMask M, typename Get>
	constexpr field_type_t<T, F, M> init_field(Get && get) noexcept(noexcept(get())) {
		if constexpr (has_field_v<M, F>) return std::forward<Get>(get)();  // call get
		else return std::monostate();
	}

	template<typename F>
	concept IsMutableFetcher =
	requires { typename std::remove_reference_t<F>::user_data_t; } &&
	IsUserData<typename std::remove_reference_t<F>::user_data_t> &&
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
		{ f.user_data() }     -> std::same_as<typename std::remove_reference_t<F>::user_data_t&>;
	} &&
	std::is_copy_constructible_v<F> &&
	std::is_move_constructible_v<F>;

	template<typename F>
	concept IsConstFetcher =
	requires { typename std::remove_reference_t<F>::user_data_t; } &&
	IsUserData<typename std::remove_reference_t<F>::user_data_t> &&
	requires(const F cf) {
		{ cf.position() }    -> std::same_as<const vec3&>;
		{ cf.velocity() }    -> std::same_as<const vec3&>;
		{ cf.force() }       -> std::same_as<const vec3&>;
		{ cf.old_position() }-> std::same_as<const vec3&>;
		{ cf.old_force() }   -> std::same_as<const vec3&>;
		{ cf.mass() }        -> std::same_as<double>;
		{ cf.state() }       -> std::same_as<ParticleState>;
		{ cf.type() }        -> std::same_as<ParticleType>;
		{ cf.id() }          -> std::same_as<ParticleID>;
		{ cf.user_data() }   -> std::same_as<const typename  std::remove_reference_t<F>::user_data_t&>;
	}  &&
	std::is_copy_constructible_v<F> &&
	std::is_move_constructible_v<F>;


	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef {

		ParticleRef(const ParticleRef& other) = default;
		ParticleRef(ParticleRef&& other) = default;

		template<class F>
		requires IsMutableFetcher<std::remove_cvref_t<F>>
		explicit ParticleRef(F && f)
			: force      ( init_field<vec3&,         		Field::force,		M>([&]() -> vec3&{ return f.force(); }) )
			, position   ( init_field<vec3&,         		Field::position,	M>([&]() -> vec3&{ return f.position(); }) )
			, velocity   ( init_field<vec3&,         		Field::velocity,	M>([&]() -> vec3&{ return f.velocity(); }) )
			, old_position(init_field<vec3&,         		Field::old_position,M>([&]() -> vec3&{ return f.old_position(); }) )
			, old_force  ( init_field<vec3&,         		Field::old_force, 	M>([&]() -> vec3&{ return f.old_force(); }) )
			, mass       ( init_field<double&,       		Field::mass,      	M>([&]() -> double&{ return f.mass(); }) )
			, state      ( init_field<ParticleState&,		Field::state,     	M>([&]() -> ParticleState&{ return f.state(); }) )
			, type       ( init_field<const ParticleType, 	Field::type,		M>([&]() -> ParticleType{ return f.type(); }) )
			, id         ( init_field<const ParticleID,   	Field::id,			M>([&]() -> ParticleID{ return f.id(); }) )
			, user_data  ( init_field<UserDataT&,			Field::user_data,	M>([&]() -> UserDataT&{ return f.user_data(); }) )
		{}

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

		RestrictedParticleRef(const RestrictedParticleRef& other) = default;
		RestrictedParticleRef(RestrictedParticleRef&& other) = default;

		template<IsMutableFetcher F>
		requires IsMutableFetcher<std::remove_cvref_t<F>>
		explicit RestrictedParticleRef(F && f) // forwarding reference
			: force        ( f.force() )
			, position     ( init_field<const vec3&,		Field::position,	M>([&]()-> const vec3&{ return f.position(); }) )
			, velocity     ( init_field<const vec3&,		Field::velocity,	M>([&]()-> const vec3&{ return f.velocity(); }) )
			, old_position ( init_field<const vec3&,		Field::old_position,M>([&]()-> const vec3&{ return f.old_position(); }) )
			, old_force    ( init_field<const vec3&,		Field::old_force,	M>([&]()-> const vec3&{ return f.old_force(); }) )
			, mass         ( init_field<const double,		Field::mass,		M>([&]()-> double { return f.mass(); }) )
			, state        ( init_field<const ParticleState,Field::state,		M>([&]()-> ParticleState { return f.state(); }) )
			, type         ( init_field<const ParticleType, Field::type,		M>([&]()-> ParticleType { return f.type(); }) )
			, id           ( init_field<const ParticleID,	Field::id,			M>([&]()-> ParticleID { return f.id(); }) )
			, user_data    ( init_field<const UserDataT&,	Field::user_data,	M>([&]()-> const UserDataT& { return f.user_data(); }) ) {
			static_assert(has_field_v<M, Field::force>,
			"ParticleRef must have Field::force to be converted to RestrictedParticleRef. "
			"Have you checked if all force fields have Field::force set?");
		}

		vec3& force;
		[[no_unique_address]] field_type_t<const vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<const vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<const double, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<const ParticleState, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<const UserDataT&, Field::user_data, M> user_data;
	};



	// Immutable reference to particle data, intended for read-only access (e.g., monitors).
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleView {

		ParticleView(const ParticleView& other) = default;
		ParticleView(ParticleView&& other) = default;

		template<class F>
		requires IsConstFetcher<std::remove_cvref_t<F>> || IsMutableFetcher<std::remove_cvref_t<F>>
		explicit ParticleView(F && f)
		   : force        ( init_field<const vec3&, 		Field::force,		M>([&]()-> const vec3&{ return f.force(); }) )
		   , position     ( init_field<const vec3&, 		Field::position,	M>([&]()-> const vec3&{ return f.position(); }) )
		   , velocity     ( init_field<const vec3&, 		Field::velocity,	M>([&]()-> const vec3&{ return f.velocity(); }) )
		   , old_position ( init_field<const vec3&, 		Field::old_position,M>([&]()-> const vec3&{ return f.old_position(); }) )
		   , old_force    ( init_field<const vec3&, 		Field::old_force,	M>([&]()-> const vec3&{ return f.old_force(); }) )
		   , mass         ( init_field<const double,		Field::mass,		M>([&]{ return f.mass(); }) )
		   , state        ( init_field<const ParticleState, Field::state,		M>([&]{ return f.state(); }) )
		   , type         ( init_field<const ParticleType,	Field::type,		M>([&]{ return f.type(); }) )
		   , id           ( init_field<const ParticleID,	Field::id,			M>([&]{ return f.id(); }) )
		   , user_data    ( init_field<const UserDataT&,	Field::user_data,	M>([&]()-> const UserDataT&{ return f.user_data(); }) )
			{}

		[[no_unique_address]] field_type_t<const vec3&, Field::force, M> force;
		[[no_unique_address]] field_type_t<const vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<const vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<const vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<const double, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<const ParticleState, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID, Field::id, M> id;
		[[no_unique_address]] field_type_t<const UserDataT&, Field::user_data, M> user_data;
	};



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



	namespace internal
	{
		// used internally in system. Holds all data of a particle
		template<IsUserData UserData>
		struct ParticleRecord {
			using user_data_t = UserData;
			ParticleRecord() = default;

			ParticleID id {};		// id of the particle.
			ParticleType type {};	// type of the particle.

			vec3 position;			// current position of the particle.
			vec3 old_position;		// previous position of the particle. Useful for applying boundary conditions
			vec3 velocity;			// current velocity of the particle.
			vec3 force;				// current force acting on the particle.
			vec3 old_force;			// previous force acting on the particle.

			ParticleState state {};	// state of the particle.
			double mass {};			// mass of the particle.

			UserData user_data; // optional user data

			bool operator==(const ParticleRecord& other) const {
				return id == other.id;
			}
		};

		template<IsUserData UserDataT>
		struct ParticleRecordFetcher {
			using user_data_t = UserDataT;
			using Record = ParticleRecord<UserDataT>;

			explicit ParticleRecordFetcher(Record& r) : record(r) {}

			vec3& position()       { return record.position; }
			vec3& velocity()       { return record.velocity; }
			vec3& force()          { return record.force; }
			vec3& old_position()   { return record.old_position; }
			vec3& old_force()      { return record.old_force; }
			double& mass()         { return record.mass; }
			ParticleState& state() { return record.state; }
			ParticleType& type()   { return record.type; }
			ParticleID& id()       { return record.id; }
			UserDataT& user_data() { return record.user_data; }
		private:
			Record& record;
		};

		template<IsUserData UserDataT>
		struct ConstParticleRecordFetcher {
			using user_data_t = UserDataT;
			using Record = ParticleRecord<UserDataT>;

			explicit ConstParticleRecordFetcher(const Record& r) : record(r) {}

			[[nodiscard]] const vec3& position() const       	{ return record.position; }
			[[nodiscard]] const vec3& velocity() const       	{ return record.velocity; }
			[[nodiscard]] const vec3& force() const          	{ return record.force; }
			[[nodiscard]] const vec3& old_position() const   	{ return record.old_position; }
			[[nodiscard]] const vec3& old_force() const      	{ return record.old_force; }
			[[nodiscard]] double mass() const         		 	{ return record.mass; }
			[[nodiscard]] ParticleState state() const 		 	{ return record.state; }
			[[nodiscard]] ParticleType type() const   		 	{ return record.type; }
			[[nodiscard]] ParticleID id() const       		 	{ return record.id; }
			[[nodiscard]] const UserDataT & user_data() const	{ return record.user_data; }

		private:
			const Record& record;
		};
	}

}


