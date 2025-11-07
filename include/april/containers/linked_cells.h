#pragma once

#include <algorithm>
#include "april/containers/contiguous_container.h"



namespace april::container {
	namespace internal {
		template <class FT, class U> class LinkedCells;
	}


	struct LinkedCells {
		template<class FT, class U> using impl = internal::LinkedCells<FT, U>;
		double cell_size_hint;

		void with_cell_size(const double cell_size) {
			cell_size_hint = cell_size;
		}
	};


	namespace internal {
		template <class FT, class U>
		class LinkedCells final : public ContiguousContainer<container::LinkedCells, FT, U> {
			using Base = ContiguousContainer<container::LinkedCells, FT, U>;
			using typename Base::ParticleRecord;
			using typename Base::ParticleID;
			using Base::force_table;
			using Base::cfg;
			using Base::domain;
			using Base::particles;
			using Base::swap_particles;
			using Base::indices;
			using Base::flags;

			enum CellWrapFlag : uint8_t {
				NO_WRAP = 0,
				WRAP_X=1,
				WRAP_Y=2,
				WRAP_Z=4,
			};

			struct CellPair {
				const uint32_t c1 = {};
				const uint32_t c2 = {};
			};

			struct WrappedCellPair {
				const uint32_t c1 = {};
				const uint32_t c2 = {};
				const CellWrapFlag force_wrap = {};
				const vec3 shift = {};
			};
		public:
			using Base::Base;

			void build(const std::vector<ParticleRecord>& input_particles) {
				this->build_storage(input_particles);
				set_cell_size();
				rebuild_cell_structure();
				compute_cell_pairs();
			}

			void register_all_particle_movements() {
				rebuild_cell_structure();
			}

			void register_particle_movement(size_t p_idx) {
				ParticleRecord & particle = particles[p_idx];

				const uint32_t dst_cell = cell_index_of(particle.position);
				// const uint32_t src_cell = cell_index_of(particle.old_position);

				// binary search to find the current cell of particle with index p_idx
				const auto it = std::ranges::upper_bound(cell_begin, p_idx);
				// Subtract one to get the cell whose range contains the index
				// std::upper_bound returns an iterator to the first element greater than value.
				const uint32_t src_cell = static_cast<uint32_t>(std::distance(cell_begin.begin(), it)) - 1;

				if (src_cell == dst_cell) return;

				++particles_per_cell[dst_cell]; // particle enters dst_cell
				--particles_per_cell[src_cell]; // particle leaves src_cell

				// destination is to the left
				if (dst_cell < src_cell) {
					for (uint32_t cell_id = src_cell; cell_id > dst_cell; --cell_id) {
						// select destination as the first (most left) position within the cell
						const size_t dst_idx = cell_begin[cell_id];

						// move particle down
						swap_particles(p_idx, dst_idx);
						++cell_begin[cell_id];

						p_idx = dst_idx;
					}
				}
				// destination is to the right
				else if (dst_cell > src_cell) {
					for (uint32_t cell_id = src_cell; cell_id < dst_cell; ++cell_id) {
						// select destination as the last (most right) position within the cell
						const size_t dst_idx = cell_begin[cell_id + 1] - 1;

						// move particle up
						swap_particles(p_idx, dst_idx);
						--cell_begin[cell_id + 1];

						p_idx = dst_idx;
					}
				}
			}

			void calculate_forces() {
				for (auto & p : particles) {
					p.old_force = p.force;
					p.force = {};
				}

				// for every cell
				for (uint32_t cid = 0; cid < cell_begin.size() - 1; cid++) {
					const uint32_t start = cell_begin[cid];
					const uint32_t end = cell_begin[cid+1];

					for (uint32_t i = start; i < end; i++) {
						for (uint32_t j = i+1; j < end; j++) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							auto pf1 = std::as_const(*this).get_fetcher_by_index(i);
							auto pf2 = std::as_const(*this).get_fetcher_by_index(j);
							vec3 force = force_table->evaluate(pf1, pf2);

							p1.force += force;
							p2.force -= force;
						}
					}
				}

				// for every neighbouring cell pair
				for (const auto & pair : neighbor_cell_pairs) {
					for (uint32_t i = cell_begin[pair.c1]; i < cell_begin[pair.c1+1]; ++i) {
						for (uint32_t j = cell_begin[pair.c2]; j < cell_begin[pair.c2+1]; ++j) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							auto pf1 = std::as_const(*this).get_fetcher_by_index(i);
							auto pf2 = std::as_const(*this).get_fetcher_by_index(j);
							vec3 force = force_table->evaluate(pf1, pf2);

							p1.force += force;
							p2.force -= force;
						}
					}
				}

				// for every wrapped cell pair
				for (const auto & pair : wrapped_cell_pairs) {
					for (uint32_t i = cell_begin[pair.c1]; i < cell_begin[pair.c1+1]; ++i) {
						for (uint32_t j = cell_begin[pair.c2]; j < cell_begin[pair.c2+1]; ++j) {
							auto & p1 = particles[i];
							auto & p2 = particles[j];

							const vec3 diff = p2.position - p1.position + pair.shift;

							auto pf1 = std::as_const(*this).get_fetcher_by_index(i);
							auto pf2 = std::as_const(*this).get_fetcher_by_index(j);
							vec3 force = force_table->evaluate(pf1, pf2, diff);

							p1.force += force;
							p2.force -= force;
						}
					}
				}
			}

			std::vector<size_t> collect_indices_in_region(const env::Box & region) {

				std::vector<uint32_t> cells = cells_in_box(region);

				std::vector<size_t> ret;

				// Reserve space for the expected average number of particles per cell
				// assuming particles are uniformly distributed
				ret.reserve(particles.size()/(cells.size()+1)); // +1 accounts for outside cell

				for (const uint32_t cid : cells) {
					const uint32_t start = cell_begin[cid];
					for (size_t i = 0; i < particles_per_cell[cid]; i++) {
						if (region.contains(particles[start + i].position) && particles[start + i].state != env::ParticleState::DEAD) {
							ret.push_back(start + i);
						}
					}
				}

				return ret;
			}


		private:
			void set_cell_size() {
				cfg.cell_size_hint = std::max(force_table->get_max_cutoff(), cfg.cell_size_hint);
				if (cfg.cell_size_hint <= 0) { // negative -> no force cutoff
					cfg.cell_size_hint = domain.extent.max();
				}

				const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.x / cfg.cell_size_hint)));
				const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.y / cfg.cell_size_hint)));
				const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.z / cfg.cell_size_hint)));

				cell_size = {domain.extent.x / num_x, domain.extent.y / num_y, domain.extent.z / num_z};
				inv_cell_size = {
					cell_size.x > 0 ? 1.0/cell_size.x : 0.0,
					cell_size.y > 0 ? 1.0/cell_size.y : 0.0,
					cell_size.z > 0 ? 1.0/cell_size.z : 0.0
				  };

				cells_per_axis = uint3{num_x, num_y, num_z};
			}

			std::vector<uint32_t> cells_in_box(const env::Box & box) {

				//  Convert world coords to cell coords (relative to domain origin)
				const vec3 min = (box.min - domain.min) * inv_cell_size;
				const vec3 max = (box.max - domain.min) * inv_cell_size;

				// clamp cell coordinates to valid ranges
				const vec3 min_clamped = {
					std::clamp(std::floor(min.x), 0.0, static_cast<double>(cells_per_axis.x - 1)),
					std::clamp(std::floor(min.y), 0.0, static_cast<double>(cells_per_axis.y - 1)),
					std::clamp(std::floor(min.z), 0.0, static_cast<double>(cells_per_axis.z - 1))
				};

				const vec3 max_clamped = {
					std::clamp(std::ceil(max.x), 0.0, static_cast<double>(cells_per_axis.x - 1)),
					std::clamp(std::ceil(max.y), 0.0, static_cast<double>(cells_per_axis.y - 1)),
					std::clamp(std::ceil(max.z), 0.0, static_cast<double>(cells_per_axis.z - 1))
				};

				const uint3 min_cell = {
					static_cast<uint32_t>(min_clamped.x),
					static_cast<uint32_t>(min_clamped.y),
					static_cast<uint32_t>(min_clamped.z)
				};

				const uint3 max_cell = {
					static_cast<uint32_t>(max_clamped.x),
					static_cast<uint32_t>(max_clamped.y),
					static_cast<uint32_t>(max_clamped.z)
				};

				std::vector<uint32_t> cells;
				for (uint32_t x = min_cell.x; x <= max_cell.x; ++x) {
					for (uint32_t y = min_cell.y; y <= max_cell.y; ++y) {
						for (uint32_t z = min_cell.z; z <= max_cell.z; ++z) {
							cells.push_back(cell_pos_to_idx(x,y,z));
						}
					}
				}

				if (!(box.min>= domain.min && box.max <= domain.max)) {
					cells.push_back(outside_cell_id);
				}

				return cells;
			}

			void rebuild_cell_structure() {
				outside_cell_id = cells_per_axis.x * cells_per_axis.y * cells_per_axis.z;
				cell_begin = std::vector<uint32_t>(cells_per_axis.x * cells_per_axis.y * cells_per_axis.z + 1);

				// how many particles in each cell
				particles_per_cell = std::vector<uint32_t>(cell_begin.size());
				for (auto & p : particles) {
					const size_t cid = cell_index_of(p.position);
					++particles_per_cell[cid];
				}

				cell_begin.assign(cell_begin.size(), 0);
				for (uint32_t c = 1; c < cell_begin.size(); ++c) {
					cell_begin[c] = cell_begin[c-1] + particles_per_cell[c-1];
				}

				// scatter particles into bins
				const size_t N = particles.size();
				std::vector<ParticleRecord> tmp_p(N);
				std::vector<uint32_t> tmp_i(N);
				std::vector<uint32_t> write_ptr = cell_begin; // copy

				for (size_t i = 0; i < N; ++i) {
					const uint32_t cid = cell_index_of(particles[i].position);
					const uint32_t dst = write_ptr[cid]++;
					tmp_p[dst] = particles[i];
					tmp_i[dst] = indices[i];
				}
				particles.swap(tmp_p);
				indices.swap(tmp_i);
			}

			void compute_cell_pairs() {
				neighbor_cell_pairs.reserve(cells_per_axis.x * cells_per_axis.y * cells_per_axis.z * 13); // rough estimate

				static const int3 displacements[13] = {
					{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
					{ 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
					{-1, 0, 1}, { 0, 1, 1}, { 0,-1, 1},
					{ 1, 1, 1}, { 1,-1, 1}, {-1, 1, 1},
					{-1,-1, 1}
				};

				auto try_wrap_cell = [&](int3 & n, vec3 & shift, int ax) -> CellWrapFlag {
					if (n[ax] == -1) {
						n[ax] = cells_per_axis[ax] - 1;
						shift[ax] = - domain.extent[ax];
					} else if (n[ax] == static_cast<int>(cells_per_axis[ax])) {
						n[ax] = 0;
						shift[ax] = domain.extent[ax];
					} else {
						return NO_WRAP;
					}
					return static_cast<CellWrapFlag>(1 << ax); // maps axis to appropriate CellWrap flag
				};

				for (const auto displacement : displacements) {
					for (unsigned int z = 0; z < cells_per_axis.z; z++) {
						for (unsigned int y = 0; y < cells_per_axis.y; y++) {
							for (unsigned int x = 0; x < cells_per_axis.x; x++) {
								const int3 base{static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)};

								int3 n = base + displacement;
								vec3 shift = {};
								int8_t wrap_flags = {};

								if (flags.periodic_x) wrap_flags |= try_wrap_cell(n, shift, 0);
								if (flags.periodic_y) wrap_flags |= try_wrap_cell(n, shift, 1);
								if (flags.periodic_z) wrap_flags |= try_wrap_cell(n, shift, 2);

								if (n.x < 0 || n.y < 0 || n.z < 0)
									continue;
								if (n.x >= static_cast<int>(cells_per_axis.x) ||
									n.y >= static_cast<int>(cells_per_axis.y) ||
									n.z >= static_cast<int>(cells_per_axis.z))
									continue;

								if (shift == vec3{}) {
									neighbor_cell_pairs.emplace_back(
										cell_pos_to_idx(x,y,z),
										cell_pos_to_idx(n.x, n.y, n.z)
									);
								} else {
									wrapped_cell_pairs.emplace_back(
										cell_pos_to_idx(x,y,z),
										cell_pos_to_idx(n.x, n.y, n.z),
										static_cast<CellWrapFlag>(wrap_flags),
										shift
									);
								}
							}
						}
					}
				}
			}

			[[nodiscard]] uint32_t cell_pos_to_idx(const uint32_t x, const uint32_t y, const uint32_t z) const noexcept{
				return  z * cells_per_axis.x * cells_per_axis.y + y * cells_per_axis.x + x;
			}

			uint32_t cell_index_of(const vec3 & position) {
				const vec3 pos = position - domain.min;
				if (pos[0] < 0 || pos[1] < 0 || pos[2] < 0) {
					return outside_cell_id;
				}

				const auto x = static_cast<uint32_t>(pos.x * inv_cell_size.x);
				const auto y = static_cast<uint32_t>(pos.y * inv_cell_size.y);
				const auto z = static_cast<uint32_t>(pos.z * inv_cell_size.z);

				if (x >= cells_per_axis.x || y >= cells_per_axis.y || z >= cells_per_axis.z) {
					return outside_cell_id;
				}

				return cell_pos_to_idx(x, y, z);
			}

			uint32_t outside_cell_id {};  // index of outside cell
			std::vector<uint32_t> cell_begin; // maps cell id to the index of the first particle within that cell
			std::vector<uint32_t> particles_per_cell; // maps cell id to the number of particles in the cell
			std::vector<CellPair> neighbor_cell_pairs;
			std::vector<WrappedCellPair> wrapped_cell_pairs;
			vec3 cell_size; // size of each cell
			vec3 inv_cell_size; // the inverse of each size component
			uint3 cells_per_axis{}; // number of cells along each axis
		};
	}
}