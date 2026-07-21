#pragma once

#include "april/containers/linked_cells/lc_core.hpp"
#include "april/containers/layout/aos.hpp"
#include "april/math/range.hpp"
#include "april/containers/batching/scalar_batch.hpp"
#include "april/containers/linked_cells/lc_batching.hpp"
#include "april/containers/linked_cells/lc_config.hpp"

namespace april::container::internal {

    template <class Config>
    class LinkedCellsAoSImpl : public LinkedCellsCore<layout::AoS<Config>> {
    	using ScalarPath =  exec::ExecutionPaths<exec::ExecutionMode::Scalar>; // AoS only supports scalar execution
    public:
    	using AsymBatch = batching::AsymmetricScalarBatch<LinkedCellsAoSImpl, ScalarPath>;
    	using SymBatch = batching::SymmetricScalarBatch<LinkedCellsAoSImpl, ScalarPath>;

    	using Base = LinkedCellsCore<layout::AoS<Config>>;

    	using Base::Base;
    	friend Base;

    	template<ParallelPolicy P, typename F>
		void for_each_interaction_batch(this auto && self, F && func) {
    		auto get_indices = [&](const size_t c, const size_t t) {
    			const size_t bin_idx = self.bin_index(c, t);
    			const size_t start = self.bin_starts[bin_idx];
    			const size_t end   = self.bin_starts[bin_idx + 1];
    			return math::Range {start, end};
    		};

    		for (const auto & phase : self.phase_schedule) {

    			self.thread_executor.template execute<P>(phase.size(), [&](size_t block_idx) {
				    thread_local LinkedCellsBatch<AsymBatch, SymBatch> batch;

					auto add_asym = [&](const math::Range & range1, const math::Range & range2) {
						AsymBatch abatch (self);
						abatch.range1 = range1;
						abatch.range2 = range2;
						batch.asym_chunks.push_back(abatch);
					};

					auto add_sym = [&](const math::Range & range) {
						SymBatch sbatch (self);
						sbatch.range = range;
						batch.sym_chunks.push_back(sbatch);
					};

					self.for_each_type_pair([&](const size_t t1, const size_t t2) {
						auto [bx, by, bz] = phase[block_idx];
						// init batch
						batch.clear();
						batch.types = {static_cast<ParticleType>(t1), static_cast<ParticleType>(t2)};

						// fill the block-batch
						self.for_each_cell_in_block(bx, by, bz, [&](size_t x, size_t y, size_t z) {
							self.process_cell_interactions(x, y, z, t1, t2, get_indices, add_sym, add_asym);
						});

						// dispatch if work exists
						if (!batch.empty()) {
							func(batch, batching::NoBatchBCP{});
						}
					});
				});

    		}

    		// handle wrapped cell pairs
    		auto process_wrapped = [&](auto&& f, const math::Range& r1, const math::Range& r2, size_t t1, size_t t2, auto&& bcp) {
    			AsymBatch wrapped_batch(self);

    			wrapped_batch.types = {static_cast<ParticleType>(t1), static_cast<ParticleType>(t2)};
    			wrapped_batch.range1 = r1;
    			wrapped_batch.range2 = r2;

    			f(wrapped_batch, bcp);
    		};

    		self.template for_each_wrapped_interaction<P>(func, get_indices, process_wrapped);
		}
    };
}


namespace april::container {
    struct LinkedCellsAoS : internal::LinkedCellsConfig{

    	template<class Config>
        using impl = internal::LinkedCellsAoSImpl<Config>;
    };
}














