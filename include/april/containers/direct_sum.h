#pragma once


#include "april/containers/container.h"
#include "april/env/particle.h"



namespace april::cont {

	namespace impl {
		template <class Env> class DirectSum;
	}

	struct DirectSum {
		template<typename  Env> using impl = impl::DirectSum<Env>;
	};

	namespace impl {
		template <class Env>
		class DirectSum final : public Container<cont::DirectSum, Env> {
			using Base = Container<cont::DirectSum, Env>;
			using typename Base::Particle;
			using typename Base::ParticleID;
			using Base::interactions;
		public:
			using Base::Base;

			void build(const std::vector<Particle> & particles);
			void calculate_forces();

			[[nodiscard]] Particle& get_particle_by_id(ParticleID) {
				throw std::runtime_error("Not implemented yet");
			}

			[[nodiscard]] ParticleID id_start() const {
				return 0;
			}

			[[nodiscard]] ParticleID id_end() const {
				return particles.size() - 1;
			}

			[[nodiscard]] Particle& get_particle_by_index(size_t index) noexcept {
				AP_ASSERT(index < particles.size(), "index must be < #particles");
				return particles[index];
			}

			[[nodiscard]] size_t index_start() const {
				return 0;
			}

			[[nodiscard]] size_t index_end() const {
				return particles.size() - 1;
			}

			[[nodiscard]] size_t particle_count() const {
				return particles.size();
			}

		private:
			std::vector<Particle> particles;
		};

		template <class Env>
		void DirectSum<Env>::build(const std::vector<Particle>& particles) {
			this->particles = std::vector(particles);
		}

		template <class Env>
		void DirectSum<Env>::calculate_forces() {
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
	}
}