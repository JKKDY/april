#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <functional>

#include "april/base/types.hpp"
#include "april/particle/defs.hpp"
#include "april/env/domain.hpp"
#include "april/particle/fields.hpp"

namespace april::container::internal {
	template <class ContainerBase>
	class LinkedCellsCore : public ContainerBase {
		using cell_index_t = uint32_t;

		enum CellWrapFlag : uint8_t {
			NO_WRAP = 0,
			WRAP_X=1,
			WRAP_Y=2,
			WRAP_Z=4,
		};

		struct CellPair {
			const cell_index_t c1 = {};
			const cell_index_t c2 = {};
		};

		struct WrappedCellPair {
			const cell_index_t c1 = {};
			const cell_index_t c2 = {};
			const CellWrapFlag force_wrap = {};
			const vec3 shift = {};
		};

	public:
		using ContainerBase::ContainerBase;
		using ContainerBase::build_storage;
		using typename ContainerBase::ParticleRecord;

		void build (this auto&& self, const std::vector<ParticleRecord>& particles) {
			self.setup_cell_grid();
			self.init_cell_order();
			self.create_neighbor_stencil();
			self.compute_wrapped_cell_pairs();
			self.build_storage(particles);
			self.pre_allocate_assignment_bins();
			self.rebuild_structure();
		}

		void rebuild_structure(this auto&& self) {
			// reset assignment vector
			for (auto& bin : self.bin_assignments) {
				bin.clear();
			}

			// repopulate assignment vector
			for (size_t i = 0; i < self.bin_sizes.size(); i++) {
				size_t start = self.bin_starts[i];
				size_t end = start + self.bin_sizes[i];
				self.template for_each_particle_view <env::Field::type | env::Field::position>(start, end,
					[&](const size_t idx, const auto & p) {
						const size_t cid = self.cell_index_from_position(p.position);
						const size_t bin = self.bin_index(cid, p.type);
						self.bin_assignments[bin].push_back(idx);
					}
				);
			}

			self.reorder_storage(self.bin_assignments);
		}

		[[nodiscard]] std::vector<size_t> collect_indices_in_region(this const auto& self, const env::Box & region) {
			std::vector<cell_index_t> cells = self.get_cells_in_region(region);
			std::vector<size_t> ret;

			// heuristic: reserve space for the expected average number of particles per cell
			const size_t est_count = cells.empty() ? 0 : (self.particle_count() * cells.size() / self.n_cells);
			ret.reserve(est_count);

			// for each cell that intersects the region: for each particle in cell perform inclusion check
			for (const size_t cid : cells) {
				const auto [start_idx, end_idx] = self.cell_index_range(cid);
				if (start_idx == end_idx) continue;

				self.template for_each_particle_view<env::Field::position | env::Field::state>(start_idx, end_idx,
					[&](const size_t i, const auto & particle) {
						if (static_cast<uint8_t>(particle.state & env::ParticleState::ALIVE) &&
							region.contains(particle.position)) {
							ret.push_back(i);
						}
					}
				);
			}

			return ret;
		}

	protected:
		size_t outside_cell_id {};
		size_t n_grid_cells {};
		size_t n_cells {}; // total cells = grid + outside
		size_t n_types {}; // types range from 0 ... n_types-1
		double global_cutoff {}; // maximum force cutoff

		vec3d cell_size; // side lengths of each cell
		vec3d inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis

		std::vector<std::vector<size_t>> bin_assignments;
		std::vector<cell_index_t> cell_ordering; // map x,y,z flat index (Nx*Ny*z+Nx*y+x) to ordering index

		// cell pair info
		std::vector<int3> neighbor_stencil;
		std::vector<WrappedCellPair> wrapped_cell_pairs;


		//------
		// SETUP
		//------
		void setup_cell_grid(this auto&& self) {
			// determine the physical cutoff (max_rc) from interactions
			double max_cutoff = 0;
			for (const auto & interaction : self.force_schema.interactions) {
				if (interaction.is_active && !interaction.used_by_types.empty() && interaction.cutoff > max_cutoff) {
					max_cutoff = interaction.cutoff;
				}
			}

			// handle edge cases (no interactions or cutoff > domain)
			// ensure at least a 2x2x2 grid (before applying cell size factor (CSF))
			if (max_cutoff <= 0 || max_cutoff > self.domain.extent.min()) {
				// TODO log warning
				max_cutoff = self.domain.extent.min() / 2.0;
			}

			double target_cell_size = self.config.get_width(max_cutoff);
			AP_ASSERT(target_cell_size > 0, "Calculated cell size must be > 0");

			// compute number of cells along each axis
			// std::floor ensures that the resulting cells are larger than or equal to target_cell_size
			const auto num_x = static_cast<cell_index_t>(std::max(1.0, std::floor(self.domain.extent.x / target_cell_size)));
			const auto num_y = static_cast<cell_index_t>(std::max(1.0, std::floor(self.domain.extent.y / target_cell_size)));
			const auto num_z = static_cast<cell_index_t>(std::max(1.0, std::floor(self.domain.extent.z / target_cell_size)));

			// calculate cell size along (stretches to fit domain exactly)
			self.cell_size = {
				self.domain.extent.x / num_x,
				self.domain.extent.y / num_y,
				self.domain.extent.z / num_z
			};

			// cache inverse (useful for fast binning: index = coord * inv_cell_size)
			self.inv_cell_size = {
				self.cell_size.x > 0 ? 1.0/self.cell_size.x : 0.0,
				self.cell_size.y > 0 ? 1.0/self.cell_size.y : 0.0,
				self.cell_size.z > 0 ? 1.0/self.cell_size.z : 0.0
			  };

			self.cells_per_axis = uint3{num_x, num_y, num_z};

			// set scalars
			self.n_types = self.force_schema.types.size();
			self.n_grid_cells = num_x * num_y * num_z;
			self.n_cells = self.n_grid_cells + 1;
			self.outside_cell_id = self.n_grid_cells;
			self.global_cutoff = max_cutoff;

			// allocate buffers
			self.bin_starts.resize(self.n_cells * self.n_types);
			self.bin_assignments.resize(self.n_cells * self.n_types);
		}

		void init_cell_order(this auto && self) {
			if (self.config.cell_ordering_fn.has_value()) {
				const std::function<std::vector<uint32_t>(uint3)> & cell_ordering_fn = self.config.cell_ordering_fn.value();
				self.cell_ordering = cell_ordering_fn(self.cells_per_axis);
			}
		}

		void create_neighbor_stencil(this auto && self) {
			// number of cells within a box of sidelengths global_cutoff
			// ceil ensures we really check every relevant cell
			const auto nx = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.x));
			const auto ny = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.y));
			const auto nz = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.z));

			const double cutoff_sq = self.global_cutoff * self.global_cutoff;

			// loop over every cell within a 2*rc sized box
			for (int z = 0; z <= nz; ++z) {  // we only need a half sphere so we can exclude cells with z<0
				for (int y = -ny; y <= ny; ++y) {
					for (int x = -nx; x <= nx; ++x) {

						// half sphere filter: only "forward" cells (use tuple ordering for filtering)
						if (std::make_tuple(z, y, x) <= std::make_tuple(0, 0, 0)) {
							continue;
						}

						// calculate distance between cell at x,y,z and center cell
						vec3d dist_vec = {
							std::abs(x) > 1 ? (std::abs(x) - 1) * self.cell_size.x : 0,
							std::abs(y) > 1 ? (std::abs(y) - 1) * self.cell_size.y : 0,
							std::abs(z) > 1 ? (std::abs(z) - 1) * self.cell_size.z : 0,
						};

						if (dist_vec.norm_squared() <= cutoff_sq) {
							self.neighbor_stencil.emplace_back(x,y,z);
						}
					}
				}
			}
		}

		void compute_wrapped_cell_pairs(this auto && self) {
			// check if cell index is out of bounds, wrap it to the other side and compute the spatial shift vector
			auto try_wrap_cell = [&](int3 & n, vec3 & shift, int ax) -> CellWrapFlag {
				const int dim_cells = static_cast<int>(self.cells_per_axis[ax]);

				if (n[ax] < 0) {
					// Wrapped "Left": Shift index by +Size
					n[ax] += dim_cells;
					shift[ax] = -self.domain.extent[ax];
				}
				else if (n[ax] >= dim_cells) {
					// Wrapped "Right": Shift index by -Size
					n[ax] -= dim_cells;
					shift[ax] = self.domain.extent[ax];
				}
				else {
					return NO_WRAP;
				}

				return static_cast<CellWrapFlag>(1 << ax); // maps axis to appropriate CellWrap flag
			};


			for (unsigned int z = 0; z < self.cells_per_axis.z; z++) {
				for (unsigned int y = 0; y < self.cells_per_axis.y; y++) {
					for (unsigned int x = 0; x < self.cells_per_axis.x; x++) {
						for (const auto displacement : self.neighbor_stencil) {
							if (displacement == int3{0,0,0}) continue;

							const int3 base {  // coordinates of cell 1
								static_cast<int>(x),
								static_cast<int>(y),
								static_cast<int>(z)
							};

							int3 n = base + displacement; // coordinates of cell 2
							vec3 shift = {};
							int8_t wrap_flags = {};

							// if periodic check if cell pair needs to be wrapped
							if (self.flags.periodic_x) wrap_flags |= try_wrap_cell(n, shift, 0);
							if (self.flags.periodic_y) wrap_flags |= try_wrap_cell(n, shift, 1);
							if (self.flags.periodic_z) wrap_flags |= try_wrap_cell(n, shift, 2);

							if (shift == vec3{}) continue;

							// if cell pair is not wrapped and cell 2 is outside of domain: skip pair
							if (n.x < 0 || n.y < 0 || n.z < 0 ||
								n.x >= static_cast<int>(self.cells_per_axis.x) ||
								n.y >= static_cast<int>(self.cells_per_axis.y) ||
								n.z >= static_cast<int>(self.cells_per_axis.z)
							) continue;

							self.wrapped_cell_pairs.emplace_back(
								self.cell_pos_to_idx(x,y,z),
								self.cell_pos_to_idx(n.x, n.y, n.z),
								static_cast<CellWrapFlag>(wrap_flags),
								shift
							);
						}
					}
				}
			}
		}

		void pre_allocate_assignment_bins() {
			const size_t num_bins = n_types * n_cells;
			bin_assignments.resize(num_bins);

			// for each bin assume somewhat uniform distribution + 50% buffer
			const size_t est_per_bin = (this->particle_count() / std::max<size_t>(1, num_bins)) * 1.5;
			for(auto& bin : bin_assignments) {
				bin.reserve(est_per_bin);
			}
		}


        // -----------------
        // LOOP ABSTRACTIONS
        // -----------------
        // Iterates over spatial blocks (cache blocking)
        template <typename Func>
        AP_FORCE_INLINE void for_each_block(Func&& fn) const {
            const auto& bdim = this->config.block_size;
            for (size_t bz = 0; bz < cells_per_axis.z; bz += bdim.z)
                for (size_t by = 0; by < cells_per_axis.y; by += bdim.y)
                    for (size_t bx = 0; bx < cells_per_axis.x; bx += bdim.x)
                        fn(bx, by, bz);
        }

		// Iterates over all unique pairs of types (T1, T2) where T2 >= T1
		template <typename Func>
		AP_FORCE_INLINE void for_each_type_pair(Func&& fn) const {
			for (size_t t1 = 0; t1 < this->n_types; ++t1)
				for (size_t t2 = t1; t2 < this->n_types; ++t2)
					fn(t1, t2);
		}

        // Iterates over the cells inside a specific block
        template <typename Func>
        AP_FORCE_INLINE void for_each_cell_in_block(size_t bx, size_t by, size_t bz, Func&& fn) const {
            const auto& bdim = this->config.block_size;

            // Calculate limits (handling edge blocks that might be smaller)
            const size_t z_end = std::min(bz + bdim.z, static_cast<size_t>(cells_per_axis.z));
            const size_t y_end = std::min(by + bdim.y, static_cast<size_t>(cells_per_axis.y));
            const size_t x_end = std::min(bx + bdim.x, static_cast<size_t>(cells_per_axis.x));

            for (size_t z = bz; z < z_end; ++z)
                for (size_t y = by; y < y_end; ++y)
                    for (size_t x = bx; x < x_end; ++x)
                        fn(x, y, z);
        }



		// ---------------------
		// BATCH ITERATOR KERNEL
		// ---------------------
		template <typename GetRange, typename AddSym, typename AddAsym>
		AP_FORCE_INLINE void process_cell_interactions(
			size_t x, size_t y, size_t z,
			size_t t1, size_t t2,
			GetRange&& get_range,
			AddSym&& add_sym,
			AddAsym&& add_asym
		) const {
			const size_t c = this->cell_pos_to_idx(x, y, z);
			auto range1 = get_range(c, t1);

			// intra-cell: process forces between particles inside the cell
			if (t1 == t2) {
				if (range1.size() > 1) add_sym(range1);
			} else {
				auto range2 = get_range(c, t2);
				if (range1.size() > 0 && range2.size() > 0) {
					add_asym(range1, range2);
				}
			}

			// If T1==T2 and Cell(T1) is empty, neighbors don't matter.
			// For mixed types (T1!=T2), we continue because we need the reverse check (Neighbor(T1) vs Cell(T2)).
			if (range1.empty() && (t1 == t2)) return;

			// inter-cell: process forces between particles of neighbouring cells
			for (auto offset : this->neighbor_stencil) {
				size_t c_n = this->get_neighbor_idx(x, y, z, offset);
				if (c_n == this->outside_cell_id) continue;

				// Interaction 1: Cell(T1) -> Neighbor(T2)
				auto range_n2 = get_range(c_n, t2);
				if (!range1.empty() && !range_n2.empty()) {
					add_asym(range1, range_n2);
				}

				// Interaction 2: Neighbor(T1) -> Cell(T2)
				// (due to the half stencil we would otherwise never iterate over this combination)
				if (t1 != t2) {
					auto range2 = get_range(c, t2);
					auto range_n1 = get_range(c_n, t1);

					if (!range2.empty() && !range_n1.empty()) {
						add_asym(range_n1, range2);
					}
				}
			}
		}

		template <typename Func, typename GetIndices, typename ProcessBatch>
		AP_FORCE_INLINE void for_each_wrapped_interaction(
			Func&& func,
			GetIndices&& get_indices,
			ProcessBatch&& process_batch
		) const {
			if (this->wrapped_cell_pairs.empty()) return;

			for (const auto& pair : this->wrapped_cell_pairs) {
				auto bcp = [&pair](const vec3& diff) { return diff + pair.shift; };

				for (size_t t1 = 0; t1 < this->n_types; ++t1) {
					auto range1 = get_indices(pair.c1, t1);
					if (range1.empty()) continue;

					for (size_t t2 = 0; t2 < this->n_types; ++t2) {
						auto range2 = get_indices(pair.c2, t2);
						if (range2.empty()) continue;

						process_batch(func, range1, range2, t1, t2, bcp);
					}
				}
			}
		}



		//----------
		// UTILITIES
		//----------
		[[nodiscard]] AP_FORCE_INLINE size_t get_neighbor_idx(size_t x, size_t y, size_t z, int3 offset) const {
			const int nx = static_cast<int>(x) + offset.x;
			const int ny = static_cast<int>(y) + offset.y;
			const int nz = static_cast<int>(z) + offset.z;

			// Note: If you implement periodic BCs via ghost cells, this check might change.
			// For now, it clamps to "Outside".
			if (nx < 0 || ny < 0 || nz < 0 ||
				nx >= static_cast<int>(this->cells_per_axis.x) ||
				ny >= static_cast<int>(this->cells_per_axis.y) ||
				nz >= static_cast<int>(this->cells_per_axis.z)) {
				return this->outside_cell_id;
				}
			return this->cell_pos_to_idx(nx, ny, nz);
		}

		// gather all cell ids whose cells have an intersection with the box region
		[[nodiscard]] std::vector<cell_index_t> get_cells_in_region(this const auto& self, const env::Box & box) {
			//  Convert world coords to cell coords (relative to domain origin)
			const vec3d min = (box.min - self.domain.min) * self.inv_cell_size;
			const vec3d max = (box.max - self.domain.min) * self.inv_cell_size;

			// clamp cell coordinates to valid ranges
			const vec3d min_clamped = {
				std::clamp(std::floor(min.x), 0.0, static_cast<double>(self.cells_per_axis.x - 1)),
				std::clamp(std::floor(min.y), 0.0, static_cast<double>(self.cells_per_axis.y - 1)),
				std::clamp(std::floor(min.z), 0.0, static_cast<double>(self.cells_per_axis.z - 1))
			};

			const vec3d max_clamped = {
				std::clamp(std::ceil(max.x), 0.0, static_cast<double>(self.cells_per_axis.x - 1)),
				std::clamp(std::ceil(max.y), 0.0, static_cast<double>(self.cells_per_axis.y - 1)),
				std::clamp(std::ceil(max.z), 0.0, static_cast<double>(self.cells_per_axis.z - 1))
			};

			// find the lowest left cell
			const uint3 min_cell = {
				static_cast<uint3::type>(min_clamped.x),
				static_cast<uint3::type>(min_clamped.y),
				static_cast<uint3::type>(min_clamped.z)
			};

			// find the highest right cell
			const uint3 max_cell = {
				static_cast<uint3::type>(max_clamped.x),
				static_cast<uint3::type>(max_clamped.y),
				static_cast<uint3::type>(max_clamped.z)
			};

			// add all cells in that region
			const auto cell_counts = max_cell - min_cell;
			std::vector<cell_index_t> cells;
			cells.reserve(cell_counts.x * cell_counts.y * cell_counts.z);
			for (size_t x = min_cell.x; x <= max_cell.x; ++x) {
				for (size_t y = min_cell.y; y <= max_cell.y; ++y) {
					for (size_t z = min_cell.z; z <= max_cell.z; ++z) {
						cells.push_back(self.cell_pos_to_idx(x,y,z));
					}
				}
			}

			if (!(box.min>= self.domain.min && box.max <= self.domain.max)) {
				cells.push_back(self.outside_cell_id);
			}

			return cells;
		}

		[[nodiscard]] size_t bin_index(const size_t cell_id, const env::ParticleType type = 0) const {
			return cell_id * n_types + static_cast<size_t>(type);
		}

		[[nodiscard]] std::pair<size_t, size_t> cell_index_range(this const auto & self, const uint32_t cid) {
			const size_t start_bin_idx = self.bin_index(cid);
			size_t start = self.bin_starts[start_bin_idx];
			size_t end = start_bin_idx + self.n_types >= self.bin_starts.size() ? self.capacity() : self.bin_starts[start_bin_idx + self.n_types];
			// size_t end = self.bin_starts[start_bin_idx + self.n_types]
			return {start, end};

		}

		[[nodiscard]] uint32_t cell_pos_to_idx(this const auto & self, const uint32_t x, const uint32_t y, const uint32_t z) noexcept{
			const uint32_t flat_idx = z * self.cells_per_axis.x * self.cells_per_axis.y + y * self.cells_per_axis.x + x;
			return self.cell_ordering.empty() ? flat_idx : self.cell_ordering[flat_idx];
		}

		uint32_t cell_index_from_position(this const auto & self, const vec3 & position) noexcept {
			const vec3 pos = position - self.domain.min;

			if (pos.x < 0 || pos.y < 0 || pos.z < 0) {
				return self.outside_cell_id;
			}

			const auto x = static_cast<uint32_t>(pos.x * self.inv_cell_size.x);
			const auto y = static_cast<uint32_t>(pos.y * self.inv_cell_size.y);
			const auto z = static_cast<uint32_t>(pos.z * self.inv_cell_size.z);

			if (x >= self.cells_per_axis.x || y >= self.cells_per_axis.y || z >= self.cells_per_axis.z) {
				return self.outside_cell_id;
			}

			return self.cell_pos_to_idx(x, y, z);
		}
	};
}
