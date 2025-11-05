#pragma once

#include "april/containers/contiguous_container.h"


namespace april::container {

	namespace internal {
		template <class FT, class U> class DirectSum;
	}

	struct DirectSum {
		template<typename FT, class U> using impl = internal::DirectSum<FT, U>;
	};

	namespace internal {
		template <class FT, class U>
		class DirectSum final : public ContiguousContainer<container::DirectSum, FT, U> {
			using Base = ContiguousContainer<container::DirectSum, FT, U>;
			using typename Base::ParticleRecord;
			using typename Base::ParticleID;
			using Base::force_table;
			using Base::particles;
			using Base::domain;
			using Base::flags;
		public:
			using Base::Base;

			void build(const std::vector<ParticleRecord> & particles) {
				this->build_storage(particles);

				// const int mode = (flags.periodic_x ? 4 : 0)
				// 	 | (flags.periodic_y ? 2 : 0)
				// 	 | (flags.periodic_z ? 1 : 0);
				// kernel = kernel_LUT[mode];
			}


			void calculate_forces()  {
				// for (auto & particle : particles) {
				// 	particle.reset_force();
				// }
				//
				// if (particles.size() < 2) return;
				// kernel(this);
			}

			std::vector<size_t> collect_indices_in_region(const env::Box & region) {
				std::vector<size_t> ret;
				// ret.reserve(static_cast<size_t>(size));

				for (size_t i = 0; i < particles.size(); i++) {
					if (region.contains(particles[i].position) && particles[i].state != env::ParticleState::DEAD) {
						ret.push_back(i);
					}
				}

				return ret;
			}

			void register_all_particle_movements() {}
			void register_particle_movement(size_t) {}

		private:
			// using KernelFn = void(*)(DirectSum*) noexcept;
			// KernelFn kernel = nullptr;
			//
			// template<bool P> static constexpr int IMIN  = P ? -1 : 0;
			// template<bool P> static constexpr int IMAX  = P ?  1 : 0;
			//
			// template<bool PX, bool PY, bool PZ>
			// static vec3 minimum_image(vec3 dr, const vec3& L) noexcept {
			// 	if constexpr (PX) dr.x -= L.x * std::round(dr.x / L.x);
			// 	if constexpr (PY) dr.y -= L.y * std::round(dr.y / L.y);
			// 	if constexpr (PZ) dr.z -= L.z * std::round(dr.z / L.z);
			// 	return dr;
			// }
			//
			// template<bool PX, bool PY, bool PZ>
			// void kernel_impl(DirectSum* self) noexcept {
			// 	const auto N = self->particles.size();
			// 	const vec3& L = self->domain.extent;
			//
			// 	for (size_t i = 0; i + 1 < N; ++i) {
			// 		auto& p1 = self->particles[i];
			// 		if (p1.state == Particle::State::DEAD) continue;
			//
			// 		for (size_t j = i + 1; j < N; ++j) {
			// 			auto& p2 = self->particles[j];
			// 			if (p2.state == Particle::State::DEAD) continue;
			//
			// 			vec3 dr = p2.position - p1.position;
			// 			dr = minimum_image<PX, PY, PZ>(dr, L);
			//
			// 			const vec3 f = self->interactions->evaluate(p1, p2, dr);
			// 			p1.force += f;
			// 			p2.force -= f; // Newton’s 3rd law
			// 		}
			// 	}
			// }
			//
			// static constexpr KernelFn kernel_LUT[8] = {
			// 	&kernel_impl<false,false,false>, // no periodicity
			// 	&kernel_impl<false,false,true >, // periodic Z
			// 	&kernel_impl<false,true ,false>, // periodic Y
			// 	&kernel_impl<false,true ,true >, // periodic YZ
			// 	&kernel_impl<true ,false,false>, // periodic X
			// 	&kernel_impl<true ,false,true >, // periodic XZ
			// 	&kernel_impl<true ,true ,false>, // periodic XY
			// 	&kernel_impl<true ,true ,true >  // periodic XYZ
			// };
			//

		};
	}
}
