#pragma once

#include "april/containers/contiguous_container.h"


namespace april::container {

	namespace internal {
		template <class Env> class DirectSum;
	}

	struct DirectSum {
		template<typename  Env> using impl = internal::DirectSum<Env>;
	};

	namespace internal {
		template <class Env>
		class DirectSum final : public ContiguousContainer<container::DirectSum, Env> {
			using Base = ContiguousContainer<container::DirectSum, Env>;
			using typename Base::Particle;
			using typename Base::ParticleID;
			using Base::interactions;
			using Base::particles;
			using Base::domain;
		public:
			using Base::Base;

			void build(const std::vector<Particle> & particles) {
				this->build_storage(particles);
			}


			void calculate_forces()  {
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

			std::vector<size_t> collect_indices_in_region(const env::Domain & region) {
				const double size = std::min(region.volume() / domain.volume(), domain.volume());

				std::vector<size_t> ret;
				ret.reserve(static_cast<size_t>(size));

				for (size_t i = 0; i < particles.size(); i++) {
					if (region.contains(particles[i].position)) {
						ret.push_back(i);
					}
				}

				return ret;
			}
		};

	}
}