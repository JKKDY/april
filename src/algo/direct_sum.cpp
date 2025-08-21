#include "april/algo/direct_sum.h"

#include "april/env/interaction.h"
#include "april/env/particle.h"

namespace april::algo::impl {

	void DirectSum::build(const std::vector<Particle> & particles) {
		this->particles = particles;
	}

	void DirectSum::calculate_forces() {
		for (auto & p : particles) {
			p.reset_force();
		}

		for (size_t i = 0; i < particles.size()-1; i++) {
			for (size_t j = i+1; j < particles.size(); j++) {
				auto & p1 = particles[i];
				auto & p2 = particles[j];

				const vec3 force = interactions.evaluate(p1, p2);

				p1.force += force;
				p2.force -= force;
			}
		}
	}
}
