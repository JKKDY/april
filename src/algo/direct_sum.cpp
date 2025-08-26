#include "april/algo/direct_sum.h"

#include "april/env/interaction.h"
#include "april/env/particle.h"

namespace april::algo::impl {

	void DirectSum::build(const std::vector<Particle> & particles) {
		this->particles = std::vector(particles);
	}

	void DirectSum::calculate_forces() {
		for (auto & particle : particles) {
			particle.reset_force();
		}

		for (size_t i = 0; i < particles.size()-1; i++) {
			for (size_t j = i+1; j < particles.size(); j++) {
				auto & p1 = particles[i];
				auto & p2 = particles[j];

				const vec3 force = interactions->evaluate(p1, p2);

				p1.force += force;
				p2.force -= force;
			}
		}
	}

	IAlgorithm::Particle& DirectSum::get_particle_by_id(ParticleID){
		throw std::runtime_error("Not implemented yet");
	}
	IAlgorithm::ParticleID DirectSum::id_start() const {
		return 0;
	}
	IAlgorithm::ParticleID DirectSum::id_end() const {
		return particles.size() - 1;
	}
	IAlgorithm::Particle& DirectSum::get_particle_by_index(const size_t index) noexcept {
		AP_ASSERT(index < particles.size(), "index must be < #particles");
		return particles[index];
	}
	size_t DirectSum::index_start() const {
		return 0;
	}
	size_t DirectSum::index_end() const {
		return particles.size() - 1;
	}

	size_t DirectSum::particle_count() const {
		return particles.size();
	}
}
