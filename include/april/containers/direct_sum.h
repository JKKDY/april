#pragma once

#include "april/containers/contiguous_container.h"


namespace april::container {

	namespace internal {
		template <class Env> class DirectSum;
	}

	struct DirectSum {
		template<typename Env> using impl = internal::DirectSum<Env>;
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
			using Base::flags;
		public:
			using Base::Base;

			void build(const std::vector<Particle> & particles) {
				this->build_storage(particles);

				const int mode = (flags.periodic_x ? 4 : 0)
					 | (flags.periodic_y ? 2 : 0)
					 | (flags.periodic_z ? 1 : 0);
				kernel = kernel_LUT[mode];
			}


			void calculate_forces()  {
				for (auto & particle : particles) {
					particle.reset_force();
				}

				if (particles.size() < 2) return;
				kernel(this);
			}

			std::vector<size_t> collect_indices_in_region(const env::Domain & region) {
				const auto box = env::Box(region);

				std::vector<size_t> ret;
				// ret.reserve(static_cast<size_t>(size));

				for (size_t i = 0; i < particles.size(); i++) {
					if (box.contains(particles[i].position) && particles[i].state != Particle::State::DEAD) {
						ret.push_back(i);
					}
				}

				return ret;
			}

			void register_all_particle_movements() {}
			void register_particle_movement(size_t) {}

		private:
			using KernelFn = void(*)(DirectSum*) noexcept;
			KernelFn kernel = nullptr;

			template<bool P> static constexpr int IMIN  = P ? -1 : 0;
			template<bool P> static constexpr int IMAX  = P ?  1 : 0;

			template<bool PX, bool PY, bool PZ>
			static void kernel_impl(DirectSum* self) noexcept {
				const auto N = self->particles.size();
				const vec3& L = self->domain.extent;

				for (size_t i = 0; i < N - 1; ++i) {
					auto& p1 = self->particles[i];
					for (size_t j = i + 1; j < N; ++j) {
						auto& p2 = self->particles[j];

						// enumerate image shifts with compile-time bounds
						for (int sx = IMIN<PX>; sx <= IMAX<PX>; ++sx)
						for (int sy = IMIN<PY>; sy <= IMAX<PY>; ++sy)
						for (int sz = IMIN<PZ>; sz <= IMAX<PZ>; ++sz) {

							const vec3 shift = vec3(sx, sy, sz) * L;
							const vec3 diff = (p2.position + shift) - p1.position;
							const vec3 f    = self->interactions->evaluate(p1, p2, diff);

							p1.force += f;
							p2.force -= f;
						}
					}
				}
			}

			static constexpr KernelFn kernel_LUT[8] = {
				&kernel_impl<false,false,false>, // no periodicity
				&kernel_impl<false,false,true >, // periodic Z
				&kernel_impl<false,true ,false>, // periodic Y
				&kernel_impl<false,true ,true >, // periodic YZ
				&kernel_impl<true ,false,false>, // periodic X
				&kernel_impl<true ,false,true >, // periodic XZ
				&kernel_impl<true ,true ,false>, // periodic XY
				&kernel_impl<true ,true ,true >  // periodic XYZ
			};


		};
	}
}

//  Total integration time: 5.15198 s
//   Total integration time: 5.11681 s
//   Total integration time: 4.93966 s
