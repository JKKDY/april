#pragma once

#include "april/containers/linked_cells/lc_core.hpp"
#include "april/containers/layout/aos.hpp"
#include "april/math/range.hpp"
#include "april/containers/batching/scalar.hpp"
#include "april/containers/linked_cells/lc_batching.hpp"
#include "april/containers/linked_cells/lc_config.hpp"

namespace april::container::internal {

    template <class Config, class U>
    class LinkedCellsAoSImpl : public LinkedCellsCore<layout::AoS<Config, U>> {
    public:
		using AsymBatch = AsymmetricScalarBatch<LinkedCellsAoSImpl>;
    	using SymBatch = SymmetricScalarBatch<LinkedCellsAoSImpl>;
		LinkedCellsBatch<AsymBatch, SymBatch> compound_batch;

    	using Base = LinkedCellsCore<layout::AoS<Config, U>>;

    	using Base::Base;
    	friend Base;

    	template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
		    const uint3 block_dim = self.config.block_size;
		    auto& batch = self.compound_batch;

    		auto add_asym_range = [&](const math::Range & range1, const math::Range & range2) {
    			AsymBatch abatch (self);
    			abatch.range1 = range1;
    			abatch.range2 = range2;
    			batch.asym_chunks.push_back(abatch);
    		};

    		auto add_sym_range = [&](const math::Range & range) {
    			SymBatch sbatch (self);
    			sbatch.range = range;
    			batch.sym_chunks.push_back(sbatch);
    		};

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

		    auto get_indices = [&](const size_t c, const size_t t) {
		        const size_t bin_idx = self.bin_index(c, t);
		        const size_t start = self.bin_starts[bin_idx];
		        const size_t end   = self.bin_starts[bin_idx + 1];
		        return math::Range {start, end};
		    };


			// LOOP ABSTRACTIONS
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


			// BATCHING KERNEL
			auto process_cell = [&](const size_t x, const size_t y, const size_t z, const size_t t1, const size_t t2) {
				const size_t c = self.cell_pos_to_idx(x, y, z);
				auto range1 = get_indices(c, t1);

				// intra-cell: process forces between particles inside the cell
				if (t1 == t2) {
					if (range1.size() > 1) add_sym_range(range1);
				} else {
					auto range2 = get_indices(c, t2);
					if (range1.size() > 0 && range2.size() > 0) {
						add_asym_range(range1, range2);
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
						add_asym_range(range1, range_n2);
					}

					// Interaction 2: Neighbor(T1) -> Cell(T2)
					// (due to the half stencil we would otherwise never iterate over this combination)
					if (t1 != t2) {
						auto range2 = get_indices(c, t2);
						auto range_n1 = get_indices(c_n, t1);

						if (!range2.empty() && !range_n1.empty()) {
							add_asym_range(range_n1, range2);
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

						AsymBatch wrapped_batch(self);
						wrapped_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
						wrapped_batch.range1 = range1;
						wrapped_batch.range2 = range2;

						func(wrapped_batch, bcp);
					}
				}
			}
		}
    };
}


namespace april::container {
    struct LinkedCellsAoS : internal::LinkedCellsConfig{
        using ConfigT = LinkedCellsAoS;

        template <class U>
        using impl = internal::LinkedCellsAoSImpl<ConfigT, U>;
    };
}