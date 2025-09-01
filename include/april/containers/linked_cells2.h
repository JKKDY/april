#pragma once

#include "april/containers/contiguous_container.h"



namespace april::cont {
	namespace impl {
		template <class Env> class LinkedCells2;
	}


	struct LinkedCells2 {
		template<typename  Env> using impl = impl::LinkedCells2<Env>;
		double cell_size_hint;
	};


	namespace impl {
		template <class Env>
		class LinkedCells2 final : public ContiguousContainer<cont::LinkedCells2, Env> {
			using Base = ContiguousContainer<cont::LinkedCells2, Env>;
			using typename Base::Particle;
			using typename Base::ParticleID;
			using Base::interactions;
			using Base::cfg;
			using Base::domain;
			using Base::particles;
			using Base::swap_particles;
			using Base::indices;
		public:
			using Base::Base;

			void build(const std::vector<Particle>& particles) {
				this->build_storage(particles);
				set_cell_size();
				build_cells();
				build_neighbour_pairs();
			}

			void calculate_forces() {
				// reset & sort particles
				build_cells();

				// for every cell
				for (uint32_t cid = 0; cid < cell_start.size() - 1; cid++) {
					const uint32_t start = cell_start[cid];
					const uint32_t end = cell_start[cid+1];

					for (uint32_t i = start; i < end; i++) {
						for (uint32_t j = i+1; j < end; j++) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							const vec3 force = interactions->evaluate(p1, p2);
							p1.force += force;
							p2.force -= force;
						}
					}
				}

				// for every cell pair
				for (auto & [c1, c2] : cell_pairs) {
					for (uint32_t i = cell_start[c1]; i < cell_start[c1+1]; i++) {
						for (uint32_t j = cell_start[c2]; j < cell_start[c2+1]; j++) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							const vec3 force = interactions->evaluate(p1, p2);
							p1.force += force;
							p2.force -= force;
						}
					}
				}
			}


		private:
			void set_cell_size() {
				cfg.cell_size_hint = std::max(interactions->get_max_cutoff(), cfg.cell_size_hint);
				if (cfg.cell_size_hint <= 0) { // negative -> no force cutoff
					cfg.cell_size_hint = domain.extent.max();
				}

				const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.x / cfg.cell_size_hint)));
				const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.y / cfg.cell_size_hint)));
				const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.z / cfg.cell_size_hint)));

				cell_extent = {domain.extent.x / num_x, domain.extent.y / num_y, domain.extent.z / num_z};
				inv_cell_extent = {
					cell_extent.x > 0 ? 1.0/cell_extent.x : 0.0,
					cell_extent.y > 0 ? 1.0/cell_extent.y : 0.0,
					cell_extent.z > 0 ? 1.0/cell_extent.z : 0.0
				  };

				num_cells = uint3{num_x, num_y, num_z};
			}

			void build_cells() {
				outside_cell = num_cells.x * num_cells.y * num_cells.z;
				cell_start = std::vector<uint32_t>(num_cells.x * num_cells.y * num_cells.z + 1);

				// how many particles in each cell
				cell_count = std::vector<uint32_t>(cell_start.size());
				for (auto & p : particles) {
					p.reset_force();
					const size_t cid = cell_of(p.position);
					++cell_count[cid];
				}

				// build cell indices
				size_t idx = 0;
				for (size_t i = 0; i < cell_start.size(); i++) {
					cell_start[i] = idx;
					idx += cell_count[i];
				}

				// scatter particles into bins
				const size_t N = particles.size();
				std::vector<Particle> tmp_p(N);
				std::vector<uint32_t> tmp_i(N);
				std::vector<uint32_t> write_ptr = cell_start; // copy

				for (size_t i = 0; i < N; ++i) {
					const uint32_t cid = cell_of(particles[i].position);
					const uint32_t dst = write_ptr[cid]++;
					tmp_p[dst] = particles[i];
					tmp_i[dst] = indices[i];
				}
				particles.swap(tmp_p);
				indices.swap(tmp_i);
			}

			void build_neighbour_pairs() {
				static const int3 displacements[13] = {
					{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
					{ 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
					{-1, 0, 1}, { 0, 1, 1}, { 0,-1, 1},
					{ 1, 1, 1}, { 1,-1, 1}, {-1, 1, 1},
					{-1,-1, 1}
				};

				for (const auto d : displacements) {
					for (unsigned int z = 0; z < num_cells.z; z++) {
						for (unsigned int y = 0; y < num_cells.y; y++) {
							for (unsigned int x = 0; x < num_cells.x; x++) {
								const int3 base{static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)};
								const int3 n = base + d;

								if (n.x < 0 || n.y < 0 || n.z < 0)
									continue;
								if (n.x >= static_cast<int>(num_cells.x) ||
									n.y >= static_cast<int>(num_cells.y) ||
									n.z >= static_cast<int>(num_cells.z))
									continue;

								cell_pairs.emplace_back(cell_pos_to_idx(x,y,z), cell_pos_to_idx(n.x, n.y, n.z));
							}
						}
					}
				}
			}

			[[nodiscard]] uint32_t cell_pos_to_idx(const uint32_t x, const uint32_t y, const uint32_t z) const noexcept{
				return  z * num_cells.x * num_cells.y + y * num_cells.x + x;
			}

			uint32_t cell_of(const vec3 & position) {
				const vec3 pos = position - domain.origin;
				if (pos[0] < 0 || pos[1] < 0 || pos[2] < 0) {
					return outside_cell;
				}

				const auto x = static_cast<uint32_t>(pos.x * inv_cell_extent.x);
				const auto y = static_cast<uint32_t>(pos.y * inv_cell_extent.y);
				const auto z = static_cast<uint32_t>(pos.z * inv_cell_extent.z);

				if (x >= num_cells.x || y >= num_cells.y || z >= num_cells.z) {
					return outside_cell;
				}

				return cell_pos_to_idx(x, y, z);
			}

			uint32_t outside_cell {};
			std::vector<uint32_t> cell_start;
			std::vector<uint32_t> cell_count;
			std::vector<std::pair<uint32_t, uint32_t>> cell_pairs;
			vec3 cell_extent;
			vec3 inv_cell_extent;
			uint3 num_cells{};
		};
	}
}