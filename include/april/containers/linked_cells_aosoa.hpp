#pragma once

#include <algorithm>
#include <ranges>

#include "linked_cells_types.hpp"
#include "april/containers/cell_orderings.hpp"
#include "april/containers/batching.hpp"
#include "april/containers/aosoa.hpp"


namespace april::container {
	namespace internal {
		template <size_t chunk_size, class U> class LinkedCellsAoSoA;
	}

	struct LinkedCellsAoSoA : internal::LinkedCellsConfig{
		template<class U> using impl = internal::LinkedCellsAoSoA<8, U>;
	};
}// namespace april::container



namespace april::container::internal {
	template <size_t chunk_size, class U>
	class LinkedCellsAoSoA : public AoSoAContainer<chunk_size, container::LinkedCellsAoSoA, U> {
		using Base = AoSoAContainer<chunk_size, container::LinkedCellsAoSoA, U>;
		friend Base;
		using typename Base::ParticleRecord;
		using Base::data;
	public:
		using Base::Base;

		//---------------
		// PUBLIC METHODS
		//---------------
		void build(this auto&& self, const std::vector<ParticleRecord> & input_particles) {
			if (self.flags.infinite_domain) {
				throw std::logic_error("infinite domain not supported on linked cells");
			}

			self.build_storage(input_particles);
			self.setup_cell_grid();
			self.init_cell_order();
			self.rebuild_structure();
			self.create_neighbor_stencil();
			self.compute_wrapped_cell_pairs();
		}

		template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
		    const uint3 block_dim = self.config.block_size;
		    auto& batch = self.compound_batch;

		    // helper with integer bounds check for getting the neighbor cell index
		    auto get_neighbor_idx = [&](const size_t x, const size_t y, const size_t z, const int3 offset) -> size_t {
		        const int nx = static_cast<int>(x) + offset.x;
		        const int ny = static_cast<int>(y) + offset.y;
		        const int nz = static_cast<int>(z) + offset.z;

		        // Check bounds (assuming non-periodic for simplicity here, or use wrapping logic)
		        if (nx < 0 || ny < 0 || nz < 0 ||
		            nx >= static_cast<int>(self.cells_per_axis.x) ||
		            ny >= static_cast<int>(self.cells_per_axis.y) ||
		            nz >= static_cast<int>(self.cells_per_axis.z)) {
		            return self.outside_cell_id;
		        }
		        return self.cell_pos_to_idx(nx, ny, nz);
		    };

			// LOOP ABSTRACTIONS
		    auto get_indices = [&](const size_t c, const size_t t) {
		        const size_t bin_idx = self.bin_index(c, t);
		        const size_t start = self.bin_start_indices[bin_idx];
		        const size_t end   = self.bin_start_indices[bin_idx + 1];
		        return std::ranges::iota_view {start, end};
		    };

			auto for_each_block = [&](auto&& fn) {
				for (size_t bz = 0; bz < self.cells_per_axis.z; bz += block_dim.z)
					for (size_t by = 0; by < self.cells_per_axis.y; by += block_dim.y)
						for (size_t bx = 0; bx < self.cells_per_axis.x; bx += block_dim.x)
							fn(bx, by, bz);
			};

			auto for_each_type_pair = [&](auto&& fn) {
				for (size_t t1 = 0; t1 < self.n_types; ++t1)
					for (size_t t2 = t1; t2 < self.n_types; ++t2)
						fn(t1, t2);
			};

			auto for_each_cell_in_block = [&](const size_t bx, const size_t by, const size_t bz, auto&& fn) {
				const size_t z_end = std::min(bz + block_dim.z, static_cast<size_t>(self.cells_per_axis.z));
				const size_t y_end = std::min(by + block_dim.y, static_cast<size_t>(self.cells_per_axis.y));
				const size_t x_end = std::min(bx + block_dim.x, static_cast<size_t>(self.cells_per_axis.x));

				for (size_t z = bz; z < z_end; ++z)
					for (size_t y = by; y < y_end; ++y)
						for (size_t x = bx; x < x_end; ++x)
							fn(x, y, z);
			};


			// KERNEL
			auto process_cell = [&](size_t x, size_t y, size_t z, size_t t1, size_t t2) {
				const size_t c = self.cell_pos_to_idx(x, y, z);
				auto range1 = get_indices(c, t1);

				// intra-cell: process forces between particles inside the cell
				if (t1 == t2) {
					if (range1.size() > 1) batch.sym_chunks.emplace_back(range1);
				} else {
					auto range2 = get_indices(c, t2);
					if (range1.size() > 0 && range2.size() > 0) {
						batch.asym_chunks.emplace_back(range1, range2);
					}
				}

				// If T1==T2 and Cell(T1) is empty, neighbors don't matter.
				// For mixed types (T1!=T2), we continue because we need the reverse check (Neighbor(T1) vs Cell(T2)).
				if (range1.empty() && (t1 == t2)) return;

				// inter-cell: process forces between particles of neighbouring cells
				for (auto offset : self.neighbor_stencil) {
					size_t c_n = get_neighbor_idx(x, y, z, offset);
					if (c_n == self.outside_cell_id) continue;

					// Interaction 1: Cell(T1) -> Neighbor(T2)
					auto range_n2 = get_indices(c_n, t2);
					if (!range1.empty() && !range_n2.empty()) {
						batch.asym_chunks.emplace_back(range1, range_n2);
					}

					// Interaction 2: Neighbor(T1) -> Cell(T2)
					// (due to the half stencil we would otherwise never iterate over this combination)
					if (t1 != t2) {
						auto range2 = get_indices(c, t2);
						auto range_n1 = get_indices(c_n, t1);

						if (!range2.empty() && !range_n1.empty()) {
							batch.asym_chunks.emplace_back(range_n1, range2);
						}
					}
				}
			};

			// EXECUTION
			for_each_block([&](size_t bx, size_t by, size_t bz) {
				for_each_type_pair([&](const size_t t1, const size_t t2) {
					// init batch
					batch.clear();
					batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};

					// fill the batch
					for_each_cell_in_block(bx, by, bz, [&](size_t x, size_t y, size_t z) {
						process_cell(x, y, z, t1, t2);
					});

					// dispatch if work exists
					if (!batch.empty()) {
						func(batch, NoBatchBCP{});
					}
				});
			});

			// handle wrapped cell pairs
			if (self.wrapped_cell_pairs.empty()) return;
			for (const auto& pair : self.wrapped_cell_pairs) {
				// define bcp (shift) function
				// TODO maybe later implement chunked batching by aggregating wrapped pairs by shift
				auto bcp = [&pair](const vec3& diff) { return diff + pair.shift; };

				for (size_t t1 = 0; t1 < self.n_types; ++t1) {
					auto range1 = get_indices(pair.c1, t1);
					if (range1.empty()) continue;

					for (size_t t2 = 0; t2 < self.n_types; ++t2) {
						auto range2 = get_indices(pair.c2, t2);
						if (range2.empty()) continue;

						using range_type = std::ranges::iota_view<size_t, size_t>;
						DirectAsymmetricBatch<range_type, range_type> wrapped_batch;
						wrapped_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
						wrapped_batch.indices1 = range1;
						wrapped_batch.indices2 = range2;

						func(wrapped_batch, bcp);
					}
				}
			}
		}


		void rebuild_structure(this auto&& self) {
			const size_t num_bins = self.bin_start_indices.size();
			std::ranges::fill(self.bin_start_indices, 0);

			// calculate the index of the first particle in each bin
			// first store the size of each bin
			for (size_t i = 0; i < self.particle_count(); i++) {
				auto p = self.template view<env::Field::type | env::Field::position>(i);
				const size_t cid = self.cell_index_from_position(p.position);

				++self.bin_start_indices[self.bin_index(cid, p.type)];
			}

			// transform "counts" -> "start indices" (index of first particle in bin)
			size_t current_sum = 0;
			for (size_t i = 0; i < num_bins; ++i) {
				const size_t count = self.bin_start_indices[i]; // read count
				self.bin_start_indices[i] = current_sum;  // write start index

				// Round UP to nearest multiple of chunk_size
				// e.g. chunk_size = 8: count=3 -> 8, count=9 -> 16.
				// const size_t capacity = (count + chunk_size - 1) & ~(chunk_size - 1);
				// current_sum += capacity;
				current_sum += count;
			}

			// self.tmp.resize(current_sum);

			// scatter particles into bins
			std::ranges::copy(self.bin_start_indices, self.write_ptr.begin());
			for (size_t i = 0; i < self.particle_count(); i++) {
				auto p = self.template view<env::Field::type | env::Field::position | env::Field::id>(i);
				const size_t cid = self.cell_index_from_position(p.position);
				const size_t dst = self.write_ptr[self.bin_index(cid, p.type)]++;

				self.write_to_tmp_storage(dst, i);
				self.id_to_index_map[p.id] = dst;
			}

			// ping-pong swap
			self.swap_tmp_storage();
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

				for (size_t i = start_idx; i < end_idx; ++i) {
					const auto & p = self.template view<env::Field::position | env::Field::state>(i);

					// TODO move particles into sentinel bucket -> avoid dead check
					if (p.state != env::ParticleState::DEAD && region.contains(p.position)) {
						ret.push_back(i);
					}
				}
			}

			return ret;
		}

	private:
		//-------------
		// DATA MEMBERS
		//-------------
		size_t outside_cell_id {};
		size_t n_grid_cells {};
		size_t n_cells {}; // total cells = grid + outside
		size_t n_types {}; // types range from 0 ... n_types-1
		double global_cutoff {}; // maximum force cutoff

		vec3d cell_size; // side lengths of each cell
		vec3d inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis

		std::vector<cell_index_t> bin_start_indices; // maps bin id to index of first particle in that bin
		std::vector<cell_index_t> cell_ordering; // map x,y,z flat index (Nx*Ny*z+Nx*y+x) to ordering index

		// used for cell rebuilding
		std::vector<size_t> write_ptr;

		// cell pair info
		std::vector<int3> neighbor_stencil;
		std::vector<WrappedCellPair> wrapped_cell_pairs;

		// batching structs
		UnifiedLCBatch compound_batch;

		ChunkedStorage<typename Base::ChunkType> tmp;


		//----------------
		// SETUP FUNCTIONS
		//----------------
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
			self.bin_start_indices.resize(self.n_cells * self.n_types + 1); // size = (total Cells * types) + sentinel
			self.write_ptr.resize(self.n_cells * self.n_types + 1);
			self.allocate_tmp_storage();
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


		//-------------
		// UTILITIES
		//-------------
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

		[[nodiscard]] std::pair<size_t, size_t> cell_index_range(const uint32_t cid) const {
			const size_t start_bin_idx = bin_index(cid);

			return {
				bin_start_indices[start_bin_idx],
				bin_start_indices[start_bin_idx + n_types]
			};
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

		 void allocate_tmp_storage() {
	       // Resize temporary storage if it's too small
	       if (tmp.n_particles < this->particle_count()) {
	          tmp.resize(this->particle_count());
	       }
	    }

	    void write_to_tmp_storage(const size_t dst_i, const size_t src_i) {
	       // 1. Locate Source in Main Data
	       //    Use the 'locate' function we just fixed (returns {chunk_idx, lane_idx})
	       const auto [src_c, src_l] = data.locate(src_i);
	       const auto& src_chunk = data.chunks[src_c];

	       // 2. Locate Destination in Tmp Data
	       const auto [dst_c, dst_l] = tmp.locate(dst_i);
	       auto& dst_chunk = tmp.chunks[dst_c];

	       // 3. Copy Data (Field by Field)
	       //    The compiler will auto-vectorize these scalar assignments since they are sequential
	       dst_chunk.pos_x[dst_l] = src_chunk.pos_x[src_l];
	       dst_chunk.pos_y[dst_l] = src_chunk.pos_y[src_l];
	       dst_chunk.pos_z[dst_l] = src_chunk.pos_z[src_l];

	       dst_chunk.vel_x[dst_l] = src_chunk.vel_x[src_l];
	       dst_chunk.vel_y[dst_l] = src_chunk.vel_y[src_l];
	       dst_chunk.vel_z[dst_l] = src_chunk.vel_z[src_l];

	       dst_chunk.frc_x[dst_l] = src_chunk.frc_x[src_l];
	       dst_chunk.frc_y[dst_l] = src_chunk.frc_y[src_l];
	       dst_chunk.frc_z[dst_l] = src_chunk.frc_z[src_l];

	       dst_chunk.old_x[dst_l] = src_chunk.old_x[src_l];
	       dst_chunk.old_y[dst_l] = src_chunk.old_y[src_l];
	       dst_chunk.old_z[dst_l] = src_chunk.old_z[src_l];

	       dst_chunk.mass[dst_l]      = src_chunk.mass[src_l];
	       dst_chunk.state[dst_l]     = src_chunk.state[src_l];
	       dst_chunk.type[dst_l]      = src_chunk.type[src_l];
	       dst_chunk.id[dst_l]        = src_chunk.id[src_l];
	       dst_chunk.user_data[dst_l] = src_chunk.user_data[src_l];
	    }

	    void swap_tmp_storage() {
	       // Efficiently swap the vectors of chunks
	       std::swap(data.chunks, tmp.chunks);
	       // Sync particle count
	       data.n_particles = tmp.n_particles;
	    }
	};
} // namespace april::container::internal


