#pragma once

#include <algorithm>
#include <ranges>

#include "linked_cells_types.hpp"
#include "april/containers/cell_orderings.hpp"
#include "april/containers/batching.hpp"
#include "april/containers/aos.hpp"
#include "april/containers/soa.hpp"



namespace april::container {
	namespace internal {
		template <class U> class LinkedCellsAoS;
		template <class U> class LinkedCellsSoA;
	}

	struct LinkedCellsAoS : internal::LinkedCellsConfig{
		template<class U> using impl = internal::LinkedCellsAoS<U>;
	};

	struct LinkedCellsSoA : internal::LinkedCellsConfig{
		template<class U> using impl = internal::LinkedCellsSoA<U>;
	};
}// namespace april::container



namespace april::container::internal {
	template <class ContainerBase>
	class LinkedCellsBase : public ContainerBase {
		friend ContainerBase;
		using typename ContainerBase::ParticleRecord;
	public:
		using ContainerBase::ContainerBase;

		//---------------
		// PUBLIC METHODS
		//---------------
		void build(this auto&& self, const std::vector<ParticleRecord> & input_particles) {
			self.build_storage(input_particles);
			self.setup_cell_grid();
			self.init_cell_order();
			self.rebuild_structure();
			self.compute_cell_pairs();

			if (self.flags.infinite_domain) {
				throw std::logic_error("infinite domain not supported on linked cells");
			}
		}

		template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
			auto get_indices = [&](const size_t cell, const env::ParticleType type) {
				const size_t bin_idx = self.bin_index(cell, type);
				const size_t start = self.bin_start_indices[bin_idx];
				const size_t end   = self.bin_start_indices[bin_idx + 1]; // +1 works because types are dense
				return std::ranges::iota_view {start, end};
			};

			// INTRA CELL
			SymmetricChunkedBatch sym_batch;
			sym_batch.chunks.reserve(self.n_grid_cells); // avoid reallocations during push back
			AsymmetricChunkedBatch asym_batch;
			asym_batch.chunks.reserve(self.n_grid_cells);
			AsymmetricChunkedBatch asym_nieghb_batch;
			asym_nieghb_batch.chunks.reserve(self.neighbor_cell_pairs.size());
			
			for (size_t t1 = 0; t1 < self.n_types; ++t1) {
				sym_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t1)};
				sym_batch.chunks.clear(); // reset size to 0 but keep capacity

				for (size_t c = 0; c < self.n_grid_cells; ++c) {
					auto range = get_indices(c, t1);
					if (range.size() < 2) continue;
					sym_batch.chunks.push_back({range});
				}

				if (!sym_batch.chunks.empty()) {
					func(sym_batch, NoBatchBCP{});
				}


				for (size_t t2 = t1 + 1; t2 < self.n_types; ++t2) {

					asym_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
					asym_batch.chunks.clear();

					for (size_t c = 0; c < self.n_grid_cells; ++c) {
						auto range1 = get_indices(c, t1);
						if (range1.empty()) continue;

						auto range2 = get_indices(c, t2);
						if (range2.empty()) continue;

						asym_batch.chunks.push_back({range1, range2});
					}

					if (!asym_batch.chunks.empty()) {
						func(asym_batch, NoBatchBCP{});
					}
				}

				// NEIGHBOR CELLS
				for (size_t t2 = 0; t2 < self.n_types; ++t2) {

					asym_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
					asym_batch.chunks.clear();

					for (const auto& pair : self.neighbor_cell_pairs) {
						auto range1 = get_indices(pair.c1, t1);
						if (range1.empty()) continue;

						auto range2 = get_indices(pair.c2, t2);
						if (range2.empty()) continue;

						asym_batch.chunks.push_back({range1, range2});
					}

					if (!asym_batch.chunks.empty()) {
						func(asym_batch, NoBatchBCP{});
					}
				}
			}


			// WRAPPED CELL PAIRS
			// (only if periodic bcp enabled)
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

						AsymmetricBatch batch;
						batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
						batch.indices1 = range1;
						batch.indices2 = range2;

						func(batch, bcp);
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
				current_sum += count;
			}

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

		vec3 cell_size; // side lengths of each cell
		vec3 inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis

		std::vector<cell_index_t> bin_start_indices; // maps bin id to index of first particle in that bin
		std::vector<cell_index_t> cell_ordering; // map x,y,z flat index (Nx*Ny*z+Nx*y+x) to ordering index

		// used for cell rebuilding
		std::vector<size_t> write_ptr;

		// cell pair info
		std::vector<int3> neighbor_stencil;
		std::vector<CellPair> neighbor_cell_pairs;
		std::vector<WrappedCellPair> wrapped_cell_pairs;


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

		void compute_cell_pairs(this auto && self) {
			// number of cells within a box of sidelengths global_cutoff
			// ceil ensures we really check every relevant cell
			const auto nx = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.x));
			const auto ny = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.y));
			const auto nz = static_cast<int>(std::ceil(self.global_cutoff * self.inv_cell_size.z));

			const double cutoff_sq = self.global_cutoff * self.global_cutoff;

			std::vector<int3> stencil;

			// loop over every cell within a 2*rc sized box
			for (int z = 0; z <= nz; ++z) {  // we only need a half sphere so we can exclude cells with z<0
				for (int y = -ny; y <= ny; ++y) {
					for (int x = -nx; x <= nx; ++x) {

						// half sphere filter: only "forward" cells (use tuple ordering for filtering)
						if (std::make_tuple(z, y, x) <= std::make_tuple(0, 0, 0)) {
							continue;
						}

						// calculate distance between cell at x,y,z and center cell
						vec3 dist_vec = {
							std::abs(x) > 1 ? (std::abs(x) - 1) * self.cell_size.x : 0,
							std::abs(y) > 1 ? (std::abs(y) - 1) * self.cell_size.y : 0,
							std::abs(z) > 1 ? (std::abs(z) - 1) * self.cell_size.z : 0,
						};

						if (dist_vec.norm_squared() <= cutoff_sq) {
							stencil.emplace_back(x,y,z);
						}
					}
				}
			}

			self.neighbor_cell_pairs.reserve(self.cells_per_axis.x * self.cells_per_axis.y * self.cells_per_axis.z * stencil.size()); // heuristic

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
						for (const auto displacement : stencil) {
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

							// if cell pair is not wrapped and cell 2 is outside of domain: skip pair
							if (n.x < 0 || n.y < 0 || n.z < 0 ||
								n.x >= static_cast<int>(self.cells_per_axis.x) ||
								n.y >= static_cast<int>(self.cells_per_axis.y) ||
								n.z >= static_cast<int>(self.cells_per_axis.z)
							) continue;

							if (shift == vec3{}) {
								self.neighbor_cell_pairs.emplace_back(
									self.cell_pos_to_idx(x,y,z),
									self.cell_pos_to_idx(n.x, n.y, n.z)
								);
							} else {
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

			self.neighbor_stencil = stencil;
		}


		//-------------
		// UTILITIES
		//-------------
		// gather all cell ids whose cells have an intersection with the box region
		[[nodiscard]] std::vector<cell_index_t> get_cells_in_region(this const auto& self, const env::Box & box) {
			//  Convert world coords to cell coords (relative to domain origin)
			const vec3 min = (box.min - self.domain.min) * self.inv_cell_size;
			const vec3 max = (box.max - self.domain.min) * self.inv_cell_size;

			// clamp cell coordinates to valid ranges
			const vec3 min_clamped = {
				std::clamp(std::floor(min.x), 0.0, static_cast<double>(self.cells_per_axis.x - 1)),
				std::clamp(std::floor(min.y), 0.0, static_cast<double>(self.cells_per_axis.y - 1)),
				std::clamp(std::floor(min.z), 0.0, static_cast<double>(self.cells_per_axis.z - 1))
			};

			const vec3 max_clamped = {
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
			uint32_t flat_idx = z * self.cells_per_axis.x * self.cells_per_axis.y + y * self.cells_per_axis.x + x;
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



	template <class U>
	class LinkedCellsAoS final : public LinkedCellsBase<AoSContainer<container::LinkedCellsAoS, U>> {
		using Base = LinkedCellsBase<AoSContainer<container::LinkedCellsAoS, U>>;
		std::vector<env::internal::ParticleRecord<U>> tmp_particles;
	public:
		using Base::particles;
		using Base::Base;

		void allocate_tmp_storage() {
			tmp_particles.resize(particles.size());
		}

		void write_to_tmp_storage(const size_t dst_idx, const size_t p_idx) {
			tmp_particles[dst_idx] = particles[p_idx];
		}

		void swap_tmp_storage() {
			std::swap(particles, tmp_particles);
		}
	};

	template <class U>
	class LinkedCellsSoA final : public LinkedCellsBase<SoAContainer<container::LinkedCellsSoA, U>> {
		using Base = LinkedCellsBase<SoAContainer<container::LinkedCellsSoA, U>>;

		// mirror the Storage type from the base class
		using Storage = typename Base::Storage;
		Storage tmp;
	public:
		using Base::Base;
		using Base::data;

		void allocate_tmp_storage() {
			if (tmp.pos_x.size() < this->particle_count()) {
				tmp.resize(this->particle_count());
			}
		}

		void write_to_tmp_storage(const size_t dst, const size_t src) {
			tmp.pos_x[dst] = data.pos_x[src]; tmp.pos_y[dst] = data.pos_y[src]; tmp.pos_z[dst] = data.pos_z[src];
			tmp.vel_x[dst] = data.vel_x[src]; tmp.vel_y[dst] = data.vel_y[src]; tmp.vel_z[dst] = data.vel_z[src];
			tmp.frc_x[dst] = data.frc_x[src]; tmp.frc_y[dst] = data.frc_y[src]; tmp.frc_z[dst] = data.frc_z[src];
			tmp.old_x[dst] = data.old_x[src]; tmp.old_y[dst] = data.old_y[src]; tmp.old_z[dst] = data.old_z[src];

			tmp.mass[dst]      = data.mass[src];
			tmp.state[dst]     = data.state[src];
			tmp.type[dst]      = data.type[src];
			tmp.id[dst]        = data.id[src];
			tmp.user_data[dst] = data.user_data[src];
		}

		void swap_tmp_storage() {
			std::swap(data.pos_x, tmp.pos_x); std::swap(data.pos_y, tmp.pos_y); std::swap(data.pos_z, tmp.pos_z);
			std::swap(data.vel_x, tmp.vel_x); std::swap(data.vel_y, tmp.vel_y); std::swap(data.vel_z, tmp.vel_z);
			std::swap(data.frc_x, tmp.frc_x); std::swap(data.frc_y, tmp.frc_y); std::swap(data.frc_z, tmp.frc_z);
			std::swap(data.old_x, tmp.old_x); std::swap(data.old_y, tmp.old_y); std::swap(data.old_z, tmp.old_z);

			std::swap(data.mass, tmp.mass);
			std::swap(data.state, tmp.state);
			std::swap(data.type, tmp.type);
			std::swap(data.id, tmp.id);
			std::swap(data.user_data, tmp.user_data);
		}
	};
} // namespace april::container::internal


