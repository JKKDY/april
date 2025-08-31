#pragma once

#include "april/env/particle.h"
#include "april/containers/contiguous_container.h"


namespace april::cont {

	namespace impl {
		template <class Env> class DirectSum;
	}

	struct DirectSum {
		template<typename  Env> using impl = impl::DirectSum<Env>;
	};

	namespace impl {
		template <class Env>
		class DirectSum final : public ContiguousContainer<cont::DirectSum, Env> {
			using Base = ContiguousContainer<cont::DirectSum, Env>;
			using typename Base::Particle;
			using typename Base::ParticleID;
			using Base::interactions;
			using Base::particles;
		public:
			using Base::Base;

			void build(const std::vector<Particle> & particles);
			void calculate_forces();
		};

		template <class Env>
		void DirectSum<Env>::build(const std::vector<Particle>& particles) {
			this->build_storage(particles);
		}

		template <class Env>
		void DirectSum<Env>::calculate_forces() {
			for (auto & particle : particles) {
				particle.reset_force();
			}

			if (particles.size() < 2) return;

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
	}
}