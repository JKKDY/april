#include "april/env/particle.h"
#include <sstream>

namespace april::env::impl {
	Particle::Particle(size_t idx, unsigned int id, const vec3& position, const vec3& velocity, double mass, unsigned int type, State state, const vec3& force, const vec3& old_force, const vec3& old_position) :
		id(id), position(position), velocity(velocity), mass(mass), type(type), state(state), force(force), old_force(old_force), old_position(old_position), index(idx) {
	}

	void Particle::update_position(const vec3& dx) {
		old_position = position;
		position += dx;
	}

	void Particle::reset_force() {
		old_force = force;
		force = vec3(0, 0, 0);
	}

	bool Particle::operator==(const Particle& other) const { 
		return index == other.index; 
	}

	std::string Particle::to_string() const {
		std::ostringstream oss;
		oss << "Particle ID: " << id << "\n"
			<< "Position: " << position[0] << position[1] << position[2] << "\n"
			<< "Velocity: " << velocity[0] << velocity[1] << velocity[2] << "\n"
			<< "Mass: " << mass << "\n"
			<< "Type: " << type << "\n"
			<< "State: " << static_cast<int>(state) << "\n";
		return oss.str();
	}
}


