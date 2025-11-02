#pragma once

#include <string>
#include <type_traits>
#include <cstdint>
#include <sstream>
#include <any>
#include <utility>

#include "april/common.h"


namespace april::env {
	enum class ParticleState : uint8_t {
		ALIVE      = 1u << 0, // Moves, exerts and experiences forces
		DEAD       = 1u << 1, // Inactive; no movement or interaction
		PASSIVE    = 1u << 2, // Moves, experiences forces but exerts none
		STATIONARY = 1u << 3, // Exerts forces but does not move or respond
		EXERTING   = ALIVE | STATIONARY, // Can exert forces on others
		MOVABLE    = ALIVE | PASSIVE,    // Can move (may or may not exert forces)
		ALL        = ~0u  // Matches all states
	};


	// Bitwise OR
	inline ParticleState operator|(ParticleState a, ParticleState b) {
		return static_cast<ParticleState>(
			static_cast<unsigned int>(a) | static_cast<unsigned int>(b)
		);
	}

	// Bitwise AND
	inline ParticleState operator&(ParticleState a, ParticleState b) {
		return static_cast<ParticleState>(
			static_cast<unsigned int>(a) & static_cast<unsigned int>(b)
		);
	}

	// Bitwise NOT
	inline ParticleState operator~(ParticleState a) {
		return static_cast<ParticleState>(
			~static_cast<unsigned int>(a)
		);
	}

	// OR-assignment
	inline ParticleState& operator|=(ParticleState& a, const ParticleState b) {
		a = a | b;
		return a;
	}

	// AND-assignment
	inline ParticleState& operator&=(ParticleState& a, const ParticleState b) {
		a = a & b;
		return a;
	}


	using ParticleType = uint16_t;
	using ParticleID = uint32_t;

	// user facing declaration with optional fields and non typed field for user data
	struct Particle {
		std::optional<ParticleID> id;			// The id of the particle.
		ParticleType type = 0;  				// The type of the particle.

		vec3 position;      					// The position of the particle.
		vec3 velocity;      					// The velocity of the particle.

		double mass{};        					// The mass of the particle.
		ParticleState state{};					// The state of the particle.

		// optional data e.g. if initializing from a simulation snapshot
		std::optional<vec3> old_position = {};		// previous position of the particle. Useful for applying boundary conditions
		std::optional<vec3> old_force = {};			// previous force acting on the particle.
		std::optional<vec3> force = {};				// current force

		std::any user_data {}; // custom user data
	};


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
	concept HasFields = requires { T::fields; };

	template<HasFields Self>
	inline constexpr FieldMask FieldOf = typename std::remove_cvref_t<Self>::fields;

	template<FieldMask M, Field F>
	inline constexpr bool has_field_v = (M & to_field_mask(F)) != 0;


	template <typename T>
	concept IsUserData =
		std::default_initializable<T> &&
		std::is_trivially_copyable_v<T> &&
		std::is_trivially_destructible_v<T> &&
		std::is_standard_layout_v<T> &&
		(!std::is_polymorphic_v<T>);

	// used to tell the environment what user data will be used
	struct NoUserData {};
	template<typename Data = NoUserData>
	struct ParticleData {
		using user_data_t = Data;
	};



	template<typename T, Field F, FieldMask M>
	using field_type_t = std::conditional_t<has_field_v<M, F>, T, std::monostate>;

	template< typename T, Field F, FieldMask M, typename Get>
	constexpr field_type_t<T, F, M> init_field(Get && get) noexcept(noexcept(get())) {
		if constexpr (has_field_v<M, F>) return std::forward<Get>(get)();  // call get
		else return std::monostate();
	}

	template<typename F, IsUserData U>
	concept IsFetcher = requires(F f, const F cf) {
		{ f.force() }         -> std::same_as<vec3&>;
		{ f.position() }      -> std::same_as<vec3&>;
		{ f.velocity() }      -> std::same_as<vec3&>;
		{ f.old_position() }  -> std::same_as<vec3&>;
		{ f.old_force() }     -> std::same_as<vec3&>;
		{ f.mass() }          -> std::same_as<double&>;
		{ f.state() }         -> std::same_as<ParticleState&>;
		{ f.type() }          -> std::same_as<ParticleType&>;
		{ f.id() }            -> std::same_as<ParticleID&>;
		{ f.user_data() }     -> std::same_as<U&>;

		{ cf.position() }	  -> std::same_as<const vec3&>;
		{ cf.velocity() }	  -> std::same_as<const vec3&>;
		{ cf.force() }		  -> std::same_as<const vec3&>;
		{ cf.old_position() } -> std::same_as<const vec3&>;
		{ cf.old_force() }	  -> std::same_as<const vec3&>;
		{ cf.mass() }		  -> std::same_as<double>;
		{ cf.state() }		  -> std::same_as<ParticleState>;
		{ cf.type() }		  -> std::same_as<ParticleType>;
		{ cf.id() }			  -> std::same_as<ParticleID>;
		{ cf.user_data() }	  -> std::same_as<const U&>;
	};

	// Reference to particle data passed to controllers and boundaries that can mutate particle data.
	template<FieldMask M, IsUserData UserDataT>
	struct ParticleRef {

		template<class Fetcher> requires IsFetcher<Fetcher, UserDataT>
		explicit ParticleRef(Fetcher& f)
			: force      ( init_field<vec3&,         		Field::force,		M>([&]{ return f.force(); }) )
			, position   ( init_field<vec3&,         		Field::position,	M>([&]{ return f.position(); }) )
			, velocity   ( init_field<vec3&,         		Field::velocity,	M>([&]{ return f.velocity(); }) )
			, old_position(init_field<vec3&,         		Field::old_position,M>([&]{ return f.old_position(); }) )
			, old_force  ( init_field<vec3&,         		Field::old_force, 	M>([&]{ return f.old_force(); }) )
			, mass       ( init_field<double&,       		Field::mass,      	M>([&]{ return f.mass(); }) )
			, state      ( init_field<ParticleState&,		Field::state,     	M>([&]{ return f.state(); }) )
			, type       ( init_field<const ParticleType&, 	Field::type,		M>([&]{ return f.type(); }) )
			, id         ( init_field<const ParticleID&,   	Field::id,			M>([&]{ return f.id(); }) )
			, user_data  ( init_field<UserDataT&,			Field::user_data,	M>([&]{ return f.user_data(); }) )
		{}

		[[no_unique_address]] field_type_t<vec3&, Field::force, M> force;
		[[no_unique_address]] field_type_t<vec3&, Field::position, M> position;
		[[no_unique_address]] field_type_t<vec3&, Field::velocity, M> velocity;
		[[no_unique_address]] field_type_t<vec3&, Field::old_position, M> old_position;
		[[no_unique_address]] field_type_t<vec3&, Field::old_force, M> old_force;
		[[no_unique_address]] field_type_t<double&, Field::mass, M> mass;
		[[no_unique_address]] field_type_t<ParticleState&, Field::state, M> state;
		[[no_unique_address]] field_type_t<const ParticleType&, Field::type, M> type;
		[[no_unique_address]] field_type_t<const ParticleID&, Field::id, M> id;
		[[no_unique_address]] field_type_t<UserDataT&, Field::user_data, M> user_data;
	};

	// Restricted reference allowing only the force field to be modified, used for fields.
	template<FieldMask M, IsUserData UserDataT>
	struct RestrictedParticleRef {
		template<class Fetcher> requires IsFetcher<Fetcher, UserDataT>
		explicit RestrictedParticleRef(Fetcher& f)
			: force        ( f.force() )
			, position     ( init_field<const vec3&,		Field::position, M>([&]{ return f.position(); }) )
			, velocity     ( init_field<const vec3&,		Field::velocity, M>([&]{ return f.velocity(); }) )
			, old_position ( init_field<const vec3&,		Field::old_position, M>([&]{ return f.old_position(); }) )
			, old_force    ( init_field<const vec3&,		Field::old_force, M>([&]{ return f.old_force(); }) )
			, mass         ( init_field<const double,		Field::mass, M>([&]{ return f.mass(); }) )
			, state        ( init_field<const ParticleState,Field::state, M>([&]{ return f.state(); }) )
			, type         ( init_field<const ParticleType, Field::type, M>([&]{ return f.type(); }) )
			, id           ( init_field<const ParticleID,	Field::id, M>([&]{ return f.id(); }) )
			, user_data    ( init_field<const UserDataT&,	Field::user_data, M>([&]{ return f.user_data(); }) )
		{}

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
		template<class Fetcher> requires IsFetcher<Fetcher, UserDataT>
		explicit ParticleView(const Fetcher& f)
		   : force        ( init_field<const vec3&, 		Field::force,		M>([&]{ return f.force(); }) )
		   , position     ( init_field<const vec3&, 		Field::position,	M>([&]{ return f.position(); }) )
		   , velocity     ( init_field<const vec3&, 		Field::velocity,	M>([&]{ return f.velocity(); }) )
		   , old_position ( init_field<const vec3&, 		Field::old_position,M>([&]{ return f.old_position(); }) )
		   , old_force    ( init_field<const vec3&, 		Field::old_force,	M>([&]{ return f.old_force(); }) )
		   , mass         ( init_field<const double,		Field::mass,		M>([&]{ return f.mass(); }) )
		   , state        ( init_field<const ParticleState, Field::state,		M>([&]{ return f.state(); }) )
		   , type         ( init_field<const ParticleType,	Field::type,		M>([&]{ return f.type(); }) )
		   , id           ( init_field<const ParticleID,	Field::id,			M>([&]{ return f.id(); }) )
		   , user_data    ( init_field<const UserDataT&,	Field::user_data,	M>([&]{ return f.user_data(); }) )
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
		template<IsUserData ParticleData>
		struct ParticleRecord {

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

			ParticleData user_data; // optional user data

			bool operator==(const ParticleRecord& other) const {
				return id == other.id;
			}
		};

		template<IsUserData UserDataT>
		struct ParticleRecordFetcher {
			using Record = ParticleRecord<UserDataT>;
			Record& record;

			explicit ParticleRecordFetcher(Record& r) : record(r) {}

			// --- Mutable accessors ---
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

			[[nodiscard]] const vec3& position() const       { return record.position; }
			[[nodiscard]] const vec3& velocity() const       { return record.velocity; }
			[[nodiscard]] const vec3& force() const          { return record.force; }
			[[nodiscard]] const vec3& old_position() const   { return record.old_position; }
			[[nodiscard]] const vec3& old_force() const      { return record.old_force; }
			[[nodiscard]] const double& mass() const         { return record.mass; }
			[[nodiscard]] const ParticleState& state() const { return record.state; }
			[[nodiscard]] const ParticleType& type() const   { return record.type; }
			[[nodiscard]] const ParticleID& id() const       { return record.id; }
			[[nodiscard]] const UserDataT& user_data() const { return record.user_data; }
		};
	}

}


