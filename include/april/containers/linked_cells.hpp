#pragma once

#include <algorithm>
#include <ranges>

#include "batch.hpp"
#include "april/containers/aos_container.hpp"



namespace april::container {
	namespace internal {
		template <class U> class LinkedCellsAoS;
		template <class U> class LinkedCellsSoA;
	}


	struct LinkedCellsAoS {
		template<class U> using impl = internal::LinkedCellsAoS<U>;
		std::optional<double> cell_size_hint;

		void with_cell_size(const double cell_size) {
			cell_size_hint = cell_size;
		}
	};

	struct LinkedCellsSoA {
		template<class U> using impl = internal::LinkedCellsSoA<U>;
		std::optional<double> cell_size_hint;

		void with_cell_size(const double cell_size) {
			cell_size_hint = cell_size;
		}
	};
}



namespace april::container::internal {
	template <class ContainerBase>
	class LinkedCellsBase : public ContainerBase {
		friend ContainerBase;
		using typename ContainerBase::ParticleRecord;


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

		struct AsymmetricChunk {
			std::ranges::iota_view<size_t, size_t> indices1;
			std::ranges::iota_view<size_t, size_t> indices2;
		};

		struct SymmetricChunk {
			std::ranges::iota_view<size_t, size_t> indices;
		};

		struct AsymmetricChunkedBatch : SerialBatch<BatchSymmetry::Asymmetric> {
			std::vector<AsymmetricChunk> chunks;
		};

		struct SymmetricChunkedBatch : SerialBatch<BatchSymmetry::Symmetric> {
			std::vector<SymmetricChunk> chunks;
		};

		struct AsymmetricBatch : SerialBatch<BatchSymmetry::Asymmetric> {
			std::ranges::iota_view<size_t, size_t> indices1;
			std::ranges::iota_view<size_t, size_t> indices2;
		};

	public:
		using ContainerBase::ContainerBase;

		void build(this auto&& self, const std::vector<ParticleRecord> & input_particles) {
			self.build_storage(input_particles);
			self.setup_cell_grid();
			self.rebuild_structure();
			self.compute_cell_pairs();

			if (self.flags.infinite_domain) {
				throw std::logic_error("infinite domain not supported on linked cells");
			}
		}


		template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
			auto get_indices = [&](const uint32_t cell, const env::ParticleType type) {
				const uint32_t bin_idx = self.bin_index(cell, type);
				const size_t start = self.bin_start_indices[bin_idx];
				const size_t end   = self.bin_start_indices[bin_idx + 1]; // +1 works because types are dense
				return std::ranges::iota_view {start, end};
			};

			// INTRA CELL
			SymmetricChunkedBatch sym_batch;
			sym_batch.chunks.reserve(self.n_grid_cells); // avoid reallocations during push back

			for (size_t t = 0; t < self.n_types; ++t) {
				sym_batch.types = {static_cast<env::ParticleType>(t), static_cast<env::ParticleType>(t)};
				sym_batch.chunks.clear(); // reset size to 0 but keep capacity

				for (uint32_t c = 0; c < self.n_grid_cells; ++c) {
					auto range = get_indices(c, t);
					if (range.size() < 2) continue;
					sym_batch.chunks.push_back({range});
				}

				if (!sym_batch.chunks.empty()) {
					func(sym_batch, NoBatchBCP{});
				}
			}

			// for every pair of types in each cell
			AsymmetricChunkedBatch asym_batch;
			asym_batch.chunks.reserve(self.n_grid_cells);

			for (size_t t1 = 0; t1 < self.n_types; ++t1) {
				for (size_t t2 = t1 + 1; t2 < self.n_types; ++t2) {

					asym_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
					asym_batch.chunks.clear();

					for (uint32_t c = 0; c < self.n_grid_cells; ++c) {
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
			}

			// NEIGHBOR CELLS
			// neighbor_cell_pairs contains pre-calculated valid pairs
			asym_batch.chunks.reserve(self.neighbor_cell_pairs.size());

			for (size_t t1 = 0; t1 < self.n_types; ++t1) {
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
			// TODO use hilbert curve sorting on the cells
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
			uint32_t current_sum = 0;
			for (size_t i = 0; i < num_bins; ++i) {
				const uint32_t count = self.bin_start_indices[i]; // read count
				self.bin_start_indices[i] = current_sum;  // write start index
				current_sum += count;
			}

			// scatter particles into bins
			std::ranges::copy(self.bin_start_indices, self.write_ptr.begin());
			for (size_t i = 0; i < self.particle_count(); i++) {
				auto p = self.template view<env::Field::type | env::Field::position | env::Field::id>(i);
				const uint32_t cid = self.cell_index_from_position(p.position);
				const uint32_t dst = self.write_ptr[self.bin_index(cid, p.type)]++;

				self.write_to_tmp_storage(dst, i);
				self.id_to_index_map[p.id] = dst;
			}

			// ping-pong swap
			self.swap_tmp_storage();
		}


		[[nodiscard]] std::vector<size_t> collect_indices_in_region(this const auto& self, const env::Box & region) {
			std::vector<uint32_t> cells = self.get_cells_in_region(region);
			std::vector<size_t> ret;

			// heuristic: reserve space for the expected average number of particles per cell
			const size_t est_count = cells.empty() ? 0 : (self.particle_count() * cells.size() / self.n_cells);
			ret.reserve(est_count);

			// for each cell that intersects the region: for each particle in cell perform inclusion check
			for (const uint32_t cid : cells) {
				const auto [start_idx, end_idx] = self.cell_index_range(cid);

				for (uint32_t i = start_idx; i < end_idx; ++i) {
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
		uint32_t outside_cell_id {};
		size_t n_grid_cells {};
		size_t n_cells {}; // total cells = grid + outside
		size_t n_types {}; // types range from 0 ... n_types-1

		vec3 cell_size; // side lengths of each cell
		vec3 inv_cell_size; // cache the inverse of each size component to avoid divisions
		uint3 cells_per_axis{}; // number of cells along each axis

		std::vector<uint32_t> bin_start_indices; // maps bin id to index of first particle in that bin

		// used for cell rebuilding
		// std::vector<ParticleRecord> tmp_particles;
		std::vector<uint32_t> write_ptr;

		// cell pair info
		std::vector<CellPair> neighbor_cell_pairs;
		std::vector<WrappedCellPair> wrapped_cell_pairs;

		
		void setup_cell_grid(this auto&& self) {
			double cell_size_hint;
			if (self.config.cell_size_hint.has_value()) {
				AP_ASSERT(config.cell_size_hint.value() > 0, "config.cell_size_hint must be greater than 0");
				cell_size_hint = self.config.cell_size_hint.value();
			} else {
				double max_cutoff = 0;
				for (const auto & interaction : self.force_schema.interactions) {
					if (interaction.is_active && !interaction.used_by_types.empty() && interaction.cutoff > max_cutoff) {
						max_cutoff = interaction.cutoff;
					}
				}

				if (max_cutoff == 0 || max_cutoff > self.domain.extent.min()) {
					max_cutoff = self.domain.extent.min() / 2;
				}

				cell_size_hint = max_cutoff;
			}

			// compute number of cells along each axis
			const auto num_x = static_cast<unsigned int>(std::max(1.0, floor(self.domain.extent.x / cell_size_hint)));
			const auto num_y = static_cast<unsigned int>(std::max(1.0, floor(self.domain.extent.y / cell_size_hint)));
			const auto num_z = static_cast<unsigned int>(std::max(1.0, floor(self.domain.extent.z / cell_size_hint)));

			// calculate cell size along each axis and cache inverse
			self.cell_size = {self.domain.extent.x / num_x, self.domain.extent.y / num_y, self.domain.extent.z / num_z};
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

			// allocate buffers
			self.bin_start_indices.resize(self.n_cells * self.n_types + 1); // size = (total Cells * types) + sentinel
			self.write_ptr.resize(self.n_cells * self.n_types + 1);
			self.allocate_tmp_storage();
		}


		void compute_cell_pairs(this auto && self) {
			self.neighbor_cell_pairs.reserve(self.cells_per_axis.x * self.cells_per_axis.y * self.cells_per_axis.z * 13); // heuristic

			static const int3 displacements[13] = {
				{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
				{ 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
				{-1, 0, 1}, { 0, 1, 1}, { 0,-1, 1},
				{ 1, 1, 1}, { 1,-1, 1}, {-1, 1, 1},
				{-1,-1, 1}
			};

			auto try_wrap_cell = [&](int3 & n, vec3 & shift, int ax) -> CellWrapFlag {
				if (n[ax] == -1) {
					n[ax] = self.cells_per_axis[ax] - 1;
					shift[ax] = - self.domain.extent[ax];
				} else if (n[ax] == static_cast<int>(self.cells_per_axis[ax])) {
					n[ax] = 0;
					shift[ax] = self.domain.extent[ax];
				} else {
					return NO_WRAP;
				}
				return static_cast<CellWrapFlag>(1 << ax); // maps axis to appropriate CellWrap flag
			};

			for (const auto displacement : displacements) {
				for (unsigned int z = 0; z < self.cells_per_axis.z; z++) {
					for (unsigned int y = 0; y < self.cells_per_axis.y; y++) {
						for (unsigned int x = 0; x < self.cells_per_axis.x; x++) {
							const int3 base{static_cast<int>(x), static_cast<int>(y), static_cast<int>(z)};

							int3 n = base + displacement;
							vec3 shift = {};
							int8_t wrap_flags = {};

							if (self.flags.periodic_x) wrap_flags |= try_wrap_cell(n, shift, 0);
							if (self.flags.periodic_y) wrap_flags |= try_wrap_cell(n, shift, 1);
							if (self.flags.periodic_z) wrap_flags |= try_wrap_cell(n, shift, 2);

							if (n.x < 0 || n.y < 0 || n.z < 0)
								continue;
							if (n.x >= static_cast<int>(self.cells_per_axis.x) ||
								n.y >= static_cast<int>(self.cells_per_axis.y) ||
								n.z >= static_cast<int>(self.cells_per_axis.z))
								continue;

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
		}


		// gather all cell ids whose cells have an intersection with the box region
		[[nodiscard]] std::vector<uint32_t> get_cells_in_region(this const auto& self, const env::Box & box) {
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
						cells.push_back(self.cell_pos_to_idx(x,y,z));
					}
				}
			}

			if (!(box.min>= self.domain.min && box.max <= self.domain.max)) {
				cells.push_back(self.outside_cell_id);
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

		uint32_t cell_index_from_position(this const auto & self, const vec3 & position) {
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

       // Temporary storage struct to mirror the SoA vectors
       struct SoABuffer {
          std::vector<double> pos_x, pos_y, pos_z;
          std::vector<double> vel_x, vel_y, vel_z;
          std::vector<double> frc_x, frc_y, frc_z;
          std::vector<double> old_x, old_y, old_z;

          std::vector<double> mass;
          std::vector<env::ParticleState> state;
          std::vector<env::ParticleType> type;
          std::vector<env::ParticleID> id;
          std::vector<U> user_data;

          void resize(size_t n) {
             pos_x.resize(n); pos_y.resize(n); pos_z.resize(n);
             vel_x.resize(n); vel_y.resize(n); vel_z.resize(n);
             frc_x.resize(n); frc_y.resize(n); frc_z.resize(n);
             old_x.resize(n); old_y.resize(n); old_z.resize(n);
             mass.resize(n); state.resize(n); type.resize(n); id.resize(n);
             user_data.resize(n);
          }
       };

       SoABuffer tmp;

    public:
       using Base::Base;

       // 1. Resize Tmp Buffer
       void allocate_tmp_storage() {
          // Only resize if we need more space
          if (tmp.pos_x.size() < this->particle_count()) {
             tmp.resize(this->particle_count());
          }
       }

       // 2. Scatter Copy (Src -> Dst)
       // We explicitly copy every field from the inherited storage to our tmp buffer
       void write_to_tmp_storage(const size_t dst, const size_t src) {
          tmp.pos_x[dst] = this->pos_x[src]; tmp.pos_y[dst] = this->pos_y[src]; tmp.pos_z[dst] = this->pos_z[src];
          tmp.vel_x[dst] = this->vel_x[src]; tmp.vel_y[dst] = this->vel_y[src]; tmp.vel_z[dst] = this->vel_z[src];
          tmp.frc_x[dst] = this->frc_x[src]; tmp.frc_y[dst] = this->frc_y[src]; tmp.frc_z[dst] = this->frc_z[src];
          tmp.old_x[dst] = this->old_x[src]; tmp.old_y[dst] = this->old_y[src]; tmp.old_z[dst] = this->old_z[src];

          tmp.mass[dst]      = this->mass[src];
          tmp.state[dst]     = this->state[src];
          tmp.type[dst]      = this->type[src];
          tmp.id[dst]        = this->id[src];
          tmp.user_data[dst] = this->user_data[src];
       }

       // 3. Commit (Swap All Vectors)
       void swap_tmp_storage() {
          std::swap(this->pos_x, tmp.pos_x); std::swap(this->pos_y, tmp.pos_y); std::swap(this->pos_z, tmp.pos_z);
          std::swap(this->vel_x, tmp.vel_x); std::swap(this->vel_y, tmp.vel_y); std::swap(this->vel_z, tmp.vel_z);
          std::swap(this->frc_x, tmp.frc_x); std::swap(this->frc_y, tmp.frc_y); std::swap(this->frc_z, tmp.frc_z);
          std::swap(this->old_x, tmp.old_x); std::swap(this->old_y, tmp.old_y); std::swap(this->old_z, tmp.old_z);

          std::swap(this->mass, tmp.mass);
          std::swap(this->state, tmp.state);
          std::swap(this->type, tmp.type);
          std::swap(this->id, tmp.id);
          std::swap(this->user_data, tmp.user_data);
       }
    };
} // namespace april::container::internal



