#pragma once

#include <algorithm>
#include <ranges>

#include "batch.h"
#include "april/containers/contiguous_container.h"



namespace april::container {
	namespace internal {
		template <class U> class LinkedCells;
	}


	struct LinkedCells {
		template<class U> using impl = internal::LinkedCells<U>;
		double cell_size_hint;

		void with_cell_size(const double cell_size) {
			cell_size_hint = cell_size;
		}
	};
}

namespace april::container::internal {
	template <class U>
	class LinkedCells final : public ContiguousContainer<container::LinkedCells, U> {
		using Base = ContiguousContainer<container::LinkedCells, U>;
		using typename Base::ParticleRecord;
		using typename Base::ParticleID;
		using Base::config;
		using Base::domain;
		using Base::particles;
		using Base::swap_particles;
		using Base::id_to_index_map;
		using Base::flags;
		using Base::init_storage;

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

		struct AsymmetricBatch : BatchBase<BatchSymmetry::asymmetric, false> {
			std::ranges::iota_view<size_t, size_t> type1_indices;
			std::ranges::iota_view<size_t, size_t> type2_indices;
		};

		struct SymmetricBatch : BatchBase<BatchSymmetry::symmetric, false> {
			std::ranges::iota_view<size_t, size_t> type_indices;
		};

	public:
		using Base::Base;

		void build(const std::vector<ParticleRecord> & input_particles) {
			init_storage(input_particles);
			setup_cell_grid();
			assign_particles_to_cells();
			compute_cell_pairs();

			if (flags.infinite_domain) {
				throw std::logic_error("infinite domain not supported on linked cells");
			}
		}


		template<typename F>
		void for_each_batch(F && func) {
			auto get_indices = [&](const uint32_t cell, const env::ParticleType type) {
				const uint32_t bin_idx = bin_index(cell, type);
				const size_t start = bin_start_indices[bin_idx];
				const size_t end   = bin_start_indices[bin_idx + 1]; // +1 works because types are dense
				return std::ranges::iota_view {start, end};
			};

			// 1. intra cell
			// for every cell, we pair every type with itself and every other type
			for (uint32_t c = 0; c < n_grid_cells; ++c) {
				for (size_t t1 = 0; t1 < n_types; ++t1) {
					auto range1 = get_indices(c, t1);
					if (range1.empty()) continue;

					SymmetricBatch batch;
					batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t1)};
					batch.type_indices = range1;
					func(batch, NoBatchBCP{});

					for (size_t t2 = t1 + 1; t2 < n_types; ++t2) {
						auto range2 = get_indices(c, t2);
						if (range2.empty()) continue;

						AsymmetricBatch abatch;
						abatch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
						abatch.type1_indices = range1;
						abatch.type2_indices = range2;

						func(abatch, NoBatchBCP{});
					}
				}
			}


			// 2. neighbor cells
			// neighbor_cell_pairs contains pre-calculated valid pairs
			for (size_t t1 = 0; t1 < n_types; ++t1) {
				for (size_t t2 = 0; t2 < n_types; ++t2) {
					for (const auto& pair : neighbor_cell_pairs) {
						auto range1 = get_indices(pair.c1, t1);
						if (range1.empty()) continue;
						auto range2 = get_indices(pair.c2, t2);
						if (range2.empty()) continue;

						AsymmetricBatch batch;
						batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
						batch.type1_indices = range1;
						batch.type2_indices = range2;

						func(batch, NoBatchBCP{});
					}
				}
			}

			// 3. wrapped cell pairs (only if periodic bcp enabled)
			// for (const auto& pair : wrapped_cell_pairs) {
			// 	// define bcp (shift) function
			// 	auto bcp = [&pair](const vec3& diff) { return diff + pair.shift; };
			//
			// 	for (size_t t1 = 0; t1 < n_types; ++t1) {
			// 		auto range1 = get_indices(pair.c1, t1);
			// 		if (range1.empty()) continue;
			//
			// 		for (size_t t2 = 0; t2 < n_types; ++t2) {
			// 			auto range2 = get_indices(pair.c2, t2);
			// 			if (range2.empty()) continue;
			//
			// 			AsymmetricBatch batch;
			// 			batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
			// 			batch.type1_indices = range1;
			// 			batch.type2_indices = range2;
			//
			// 			func(batch, bcp);
			// 		}
			// 	}
			// }
		}


		void register_all_particle_movements() {
			assign_particles_to_cells();
		}

		void register_particle_movement(size_t) {
			assign_particles_to_cells();
		}

		std::vector<size_t> collect_indices_in_region(const env::Box & region) {
			std::vector<uint32_t> cells = get_cells_in_region(region);
			std::vector<size_t> ret;

			// heuristic: reserve space for the expected average number of particles per cell
			const size_t est_count = cells.empty() ? 0 : (particles.size() * cells.size() / n_cells);
			ret.reserve(est_count);

			// for each cell that intersects the region: for each particle in cell perform inclusion check
			for (const uint32_t cid : cells) {
				const auto [start_idx, end_idx] = cell_index_range(cid);

				for (uint32_t i = start_idx; i < end_idx; ++i) {
					const auto& p = particles[i];

					// TODO move particles into sentinel bucket -> avoid dead check
					if (p.state != env::ParticleState::DEAD && region.contains(p.position)) {
						ret.push_back(i);
					}
				}
			}

			return ret;
		}

	private:
		uint32_t outside_cell_id {};
		size_t n_grid_cells {};
		size_t n_cells {}; // total cells = grid + outside
		size_t n_types {}; // types range from 0 ... n_types-1

		vec3 cell_size; // side lengths of each cell
		vec3 inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis

		std::vector<uint32_t> bin_start_indices; // maps bin id to index of first particle in that bin

		// used for cell rebuilding
		std::vector<ParticleRecord> tmp_particles;
		std::vector<uint32_t> write_ptr;

		// cell pair info
		std::vector<CellPair> neighbor_cell_pairs;
		std::vector<WrappedCellPair> wrapped_cell_pairs;

		
		void setup_cell_grid() {
			// TODO add automatic cell size option

			// compute number of cells along each axis
			const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.x / config.cell_size_hint)));
			const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.y / config.cell_size_hint)));
			const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(domain.extent.z / config.cell_size_hint)));

			// calculate cell size along each axis and cache inverse
			cell_size = {domain.extent.x / num_x, domain.extent.y / num_y, domain.extent.z / num_z};
			inv_cell_size = {
				cell_size.x > 0 ? 1.0/cell_size.x : 0.0,
				cell_size.y > 0 ? 1.0/cell_size.y : 0.0,
				cell_size.z > 0 ? 1.0/cell_size.z : 0.0
			  };

			cells_per_axis = uint3{num_x, num_y, num_z};

			// compute number of types
			// TODO this info should be injected by ContainerHints
			std::unordered_set<env::ParticleType> types;
			for (const auto & p : particles) {
				types.insert(p.type);
			}

			// set scalars
			n_types = types.size();
			n_grid_cells = num_x * num_y * num_z;
			n_cells = n_grid_cells + 1;
			outside_cell_id = n_grid_cells;

			// allocate buffers
			bin_start_indices.resize(n_cells * n_types + 1); // size = (total Cells * types) + sentinel
			write_ptr.resize(n_cells * n_types);
			tmp_particles.resize(particles.size());
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

		void assign_particles_to_cells() {
			const size_t num_bins = bin_start_indices.size();
			std::ranges::fill(bin_start_indices, 0);

			// calculate the index of the first particle in each bin
			// first store the size of each bin
			for (const auto & p : particles) {
				const env::ParticleType type_idx = p.type;
				const size_t cid = cell_index_from_position(p.position);

				++bin_start_indices[bin_index(cid, type_idx)];
			}

			// transform "counts" -> "start indices" (index of first particle in bin)
			uint32_t current_sum = 0;
			for (size_t i = 0; i < num_bins; ++i) {
				const uint32_t count = bin_start_indices[i]; // read count
				bin_start_indices[i] = current_sum;  // write start index
				current_sum += count;
			}

			// scatter particles into bins
			std::ranges::copy(bin_start_indices, write_ptr.begin());
			for (const auto & p : particles) {
				const uint32_t cid = cell_index_from_position(p.position);
				const uint32_t dst = write_ptr[bin_index(cid, p.type)]++;

				tmp_particles[dst] = p;
				id_to_index_map[p.id] = dst;
			}

			// ping-pong swap
			std::swap(particles, tmp_particles);
		}


		// gather all cell ids whose cells have an intersection with the box region
		std::vector<uint32_t> get_cells_in_region(const env::Box & box) {
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

			// find the lowest left cell
			const uint3 min_cell = {
				static_cast<uint32_t>(min_clamped.x),
				static_cast<uint32_t>(min_clamped.y),
				static_cast<uint32_t>(min_clamped.z)
			};

			// find the highest right cell
			const uint3 max_cell = {
				static_cast<uint32_t>(max_clamped.x),
				static_cast<uint32_t>(max_clamped.y),
				static_cast<uint32_t>(max_clamped.z)
			};

			// add all cells in that region
			const auto cell_counts = max_cell - min_cell;
			std::vector<uint32_t> cells;
			cells.reserve(cell_counts.x * cell_counts.y * cell_counts.z);
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


		// ----- Utilities -----
		[[nodiscard]] size_t bin_index(const size_t cell_id, const env::ParticleType type = 0) const {
			return cell_id * n_types + static_cast<size_t>(type);
		}

		[[nodiscard]] std::pair<uint32_t, uint32_t> cell_index_range(const uint32_t cid) const {
			const size_t start_bin_idx = bin_index(cid);

			return {
				bin_start_indices[start_bin_idx],
				bin_start_indices[start_bin_idx + n_types]
			};
		}

		[[nodiscard]] uint32_t cell_pos_to_idx(const uint32_t x, const uint32_t y, const uint32_t z) const noexcept{
			return  z * cells_per_axis.x * cells_per_axis.y + y * cells_per_axis.x + x;
		}

		uint32_t cell_index_from_position(const vec3 & position) {
			const vec3 pos = position - domain.min;
			if (pos.x < 0 || pos.y < 0 || pos.z < 0) {
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
	};
}


