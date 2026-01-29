#pragma once

#include "april/containers/linked_cells/lc_core.hpp"
#include "april/containers/layout/aosoa.hpp"
#include "april/math/range.hpp"
#include "april/containers/batching/chunked.hpp"
#include "april/containers/linked_cells/lc_batching.hpp"
#include "april/containers/linked_cells/lc_config.hpp"

namespace april::container::internal {

    template <class Config, class U, size_t ChunkSize>
    class LinkedCellsAoSoAImpl : public LinkedCellsCore<layout::AoSoA<Config, U, ChunkSize>> {
    public:
    	using Base = LinkedCellsCore<layout::AoSoA<Config, U>>;
		using AsymBatch = AsymmetricChunkedBatch<LinkedCellsAoSoAImpl, typename Base::ChunkT>;
    	using SymBatch = SymmetricChunkedBatch<LinkedCellsAoSoAImpl, typename Base::ChunkT>;

    	using Base::Base;
    	friend Base;

    	template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
		    // const uint3 block_dim = self.config.block_size;

    		struct BinRange {
    			math::Range range_chunks;
    			size_t tail{};
    			size_t n_particles{};

    			[[nodiscard]] size_t size() const {
    				return n_particles;
    			}

    			[[nodiscard]] size_t empty() const {
    				return n_particles == 0;
    			}
    		};

		    auto get_indices = [&](const size_t c, const size_t t) {
		        const size_t bin_idx = self.bin_index(c, t);
		    	const size_t start = self.bin_starts[bin_idx];
		    	const size_t end = self.bin_starts[bin_idx+1];
		    	const size_t size = self.bin_sizes[bin_idx];

		    	const math::Range chunks {start >> self.chunk_shift, end >> self.chunk_shift};
		    	const size_t tail = size & self.chunk_mask;
		        return BinRange {chunks, tail, size};
		    };

    		auto add_asym = [&](const BinRange & range1, const BinRange & range2) {
    			AsymBatch ab (self, self.ptr_chunks);
    			ab.range1_chunks = range1.range_chunks;
    			ab.range2_chunks = range2.range_chunks;
    			ab.range1_tail = range1.tail;
    			ab.range2_tail = range2.tail;
    			self.batch.asym_chunks.push_back(ab);
    		};

    		auto add_sym = [&](const BinRange & range) {
    			SymBatch sb (self, self.ptr_chunks);
    			sb.range_chunks = range.range_chunks;
    			sb.range_tail = range.tail;
    			self.batch.sym_chunks.push_back(sb);
    		};


			// EXECUTION
			self.for_each_block([&](size_t bx, size_t by, size_t bz) {
				self.for_each_type_pair([&](const size_t t1, const size_t t2) {
					// init batch
					self.batch.clear();
					self.batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};

					// fill the batch
					self.for_each_cell_in_block(bx, by, bz, [&](size_t x, size_t y, size_t z) {
						self.process_cell_interactions(x, y, z, t1, t2,
							get_indices, add_sym, add_asym);
					});

					// dispatch if work exists
					if (!self.batch.empty()) {
						func(self.batch, NoBatchBCP{});
					}
				});
			});

			// handle wrapped cell pairs
    		auto process_wrapped = [&](auto&& f, const BinRange& r1, const BinRange& r2, size_t t1, size_t t2, auto&& bcp) {
    			AsymBatch wrapped_batch(self, self.ptr_chunks);

    			wrapped_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
    			wrapped_batch.range1_chunks = r1.range_chunks;
    			wrapped_batch.range1_tail   = r1.tail;
    			wrapped_batch.range2_chunks = r2.range_chunks;
    			wrapped_batch.range2_tail   = r2.tail;

    			f(wrapped_batch, bcp);
    		};

    		self.for_each_wrapped_interaction(func, get_indices, process_wrapped);
		}

    private:
    	LinkedCellsBatch<AsymBatch, SymBatch> batch;
    };
}


namespace april::container {

    struct LinkedCellsAoSoA : internal::LinkedCellsConfig{
        using ConfigT = LinkedCellsAoSoA;

        template <class U>
        using impl = internal::LinkedCellsAoSoAImpl<ConfigT, U, 8>;
    };

	template<size_t ChunkSize>
	struct LinkedCellsAoSoA_withChunkSize : internal::LinkedCellsConfig{
		using ConfigT = LinkedCellsAoSoA_withChunkSize;

		template <class U>
		using impl = internal::LinkedCellsAoSoAImpl<ConfigT, U, ChunkSize>;
	};
}