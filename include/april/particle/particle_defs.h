#pragma once


namespace april::env {

	enum class ParticleState : uint8_t {
		ALIVE      = 1u << 0, // Moves, exerts and experiences forces
		DEAD       = 1u << 1, // Inactive; no movement or interaction
		PASSIVE    = 1u << 2, // Moves, experiences forces but exerts none
		STATIONARY = 1u << 3, // Exerts forces but does not move or respond
		EXERTING   = ALIVE | STATIONARY, // Can exert forces on others
		MOVABLE    = ALIVE | PASSIVE,    // Can move (may or may not exert forces)
		ALL        = 0b11111111  // Matches all states
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


	template <typename T>
		concept IsUserData =
			std::default_initializable<T> &&
			std::is_trivially_copyable_v<T> &&
			std::is_trivially_destructible_v<T> &&
			std::is_standard_layout_v<T> &&
			(!std::is_polymorphic_v<T>);


	struct NoUserData {};

	// used to tell the environment what user data will be used
	template<typename Data = NoUserData>
	struct ParticleData {
		using user_data_t = Data;
	};

	template<typename Data = NoUserData>
	inline constexpr ParticleData<Data> particle_data {};

}
