#include "april/env/particle.h"
#include <sstream>



namespace april::env::internal {
	Particle::Particle(const unsigned int id, const vec3& position, const vec3& velocity,
					   const double mass, const unsigned int type, const State state, const vec3& force,
					   const vec3& old_force, const vec3& old_position) :
		position(position)
		, old_position(old_position)
		, velocity(velocity)
		, force(force)
		, old_force(old_force)
		, state(state)
		, mass(mass)
		, type(type)
		, id(id) {}


	void Particle::update_position(const vec3& dx) noexcept {
		old_position = position;
		position += dx;
	}

	void Particle::update_velocity(const vec3& dv) noexcept {
		velocity += dv;
	}

	void Particle::update_force(const vec3& df) noexcept {
		force += df;
	}

	void Particle::reset_force() noexcept {
		old_force = force;
		force = vec3(0, 0, 0);
	}

	bool Particle::operator==(const Particle& other) const {
		return id == other.id;
	}

	std::string Particle::to_string() const {
		std::ostringstream oss;
		oss << "Particle ID: " << id << "\n"
			<< "Position: " << position[0] << " " << position[1] << " " << position[2] << "\n"
			<< "Velocity: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n"
			<< "Force: " << force[0] << " " << force[1] << " " << force[2] << "\n"
			<< "Mass: " << mass << "\n"
			<< "Type: " << type << "\n"
			<< "State: " << static_cast<int>(state) << "\n";
		return oss.str();
	}
}

namespace april::env {

	ParticleRef::ParticleRef(internal::Particle& p)
		: position(p.position)
		, old_position(p.old_position)
		, velocity(p.velocity)
		, force(p.force)
		, old_force(p.old_force)
		, state(p.state)
		, mass(p.mass)
		, type(p.type)
		, id(p.id) {}

	void ParticleRef::update_position(const vec3& dx) const noexcept {
		old_position = position;
		position += dx;
	}

	void ParticleRef::update_velocity(const vec3& dv) const noexcept {
		velocity += dv;
	}

	void ParticleRef::update_force(const vec3& df) const noexcept {
		force += df;
	}

	void ParticleRef::reset_force() const noexcept {
		old_force = force;
		force = vec3(0, 0, 0);
	}

	bool ParticleRef::operator==(const internal::Particle& other) const {
		return id == other.id;
	}

	std::string ParticleRef::to_string() const {
		std::ostringstream oss;
		oss << "Particle ID: " << id << "\n"
			<< "Position: " << position[0] << " " << position[1] << " " << position[2] << "\n"
			<< "Velocity: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n"
			<< "Force: " << force[0] << " " << force[1] << " " << force[2] << "\n"
			<< "Mass: " << mass << "\n"
			<< "Type: " << type << "\n"
			<< "State: " << static_cast<int>(state) << "\n";
		return oss.str();
	}


	ParticleView::ParticleView(const internal::Particle& p):
	position(p.position)
	, old_position(p.old_position)
	, velocity(p.velocity)
	, force(p.force)
	, old_force(p.old_force)
	, state(p.state)
	, mass(p.mass)
	, type(p.type)
	, id(p.id) {}

	bool ParticleView::operator==(const internal::Particle& other) const {
		return id == other.id;
	}

	std::string ParticleView::to_string() const {
		std::ostringstream oss;
		oss << "Particle ID: " << id << "\n"
			<< "Position: " << position[0] << " " << position[1] << " " << position[2] << "\n"
			<< "Velocity: " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n"
			<< "Force: " << force[0] << " " << force[1] << " " << force[2] << "\n"
			<< "Mass: " << mass << "\n"
			<< "Type: " << type << "\n"
			<< "State: " << static_cast<int>(state) << "\n";
		return oss.str();
	}
}
