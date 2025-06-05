#include "april/env/particle.h"
#include <sstream>


namespace april::env {
	// Particle::Particle(const vec3& position, const vec3& velocity, const double mass, const ParticleType type, const ParticleID id,
	// 						   const ParticleState state):
	// 	id(id), type(type), position((position)), velocity(velocity), mass(mass), state(state) {}
}

namespace april::env::impl {
	Particle::Particle(const size_t index, const unsigned int id, const vec3& position, const vec3& velocity,
	                   const double mass, const unsigned int type, const State state, const vec3& force,
	                   const vec3& old_force, const vec3& old_position) :
		position(position), old_position(old_position), velocity(velocity), force(force), old_force(old_force),
		state(state), mass(mass), type(type), id(id), index(index) {}

	bool Particle::operator==(const Particle& other) const {
		return index == other.index;
	}

	std::string Particle::to_string() const {
		std::ostringstream oss;
		oss << "Particle ID: " << id << "\n"
			<< "Position: " << position[0] << " " << position[1] << " " << position[2] << "\n"
			<< "Velocity: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n"
			<< "Mass: " << mass << "\n"
			<< "Type: " << type << "\n"
			<< "State: " << static_cast<int>(state) << "\n";
		return oss.str();
	}
}


