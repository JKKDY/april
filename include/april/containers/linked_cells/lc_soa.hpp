#pragma once

#include "april/containers/linked_cells/lc_core.hpp"
#include "april/containers/layout/soa.hpp"
#include "april/math/range.hpp"
#include "april/containers/batching/scalar.hpp"
#include "april/containers/linked_cells/lc_batching.hpp"
#include "april/containers/linked_cells/lc_config.hpp"

namespace april::container::internal {

    template <class Config, class U>
    class LinkedCellsSoAImpl : public LinkedCellsCore<layout::SoA<Config, U>> {
    public:
		using AsymBatch = AsymmetricScalarBatch<LinkedCellsSoAImpl>;
    	using SymBatch = SymmetricScalarBatch<LinkedCellsSoAImpl>;

    	using Base = LinkedCellsCore<layout::SoA<Config, U>>;

    	using Base::Base;
    	friend Base;

    	template<typename F>
		void for_each_interaction_batch(this auto && self, F && func) {

    		auto add_asym = [&](const math::Range & range1, const math::Range & range2) {
    			AsymBatch abatch (self);
    			abatch.range1 = range1;
    			abatch.range2 = range2;
    			self.batch.asym_chunks.push_back(abatch);
    		};

    		auto add_sym = [&](const math::Range & range) {
    			SymBatch sbatch (self);
    			sbatch.range = range;
    			self.batch.sym_chunks.push_back(sbatch);
    		};

		    auto get_indices = [&](const size_t c, const size_t t) {
		        const size_t bin_idx = self.bin_index(c, t);
		        const size_t start = self.bin_starts[bin_idx];
		        const size_t end   = self.bin_starts[bin_idx + 1];
		        return math::Range {start, end};
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
    		auto process_wrapped = [&](auto&& f, const math::Range& r1, const math::Range& r2, size_t t1, size_t t2, auto&& bcp) {
    			AsymBatch wrapped_batch(self);

    			wrapped_batch.types = {static_cast<env::ParticleType>(t1), static_cast<env::ParticleType>(t2)};
    			wrapped_batch.range1 = r1;
    			wrapped_batch.range2 = r2;

    			f(wrapped_batch, bcp);
    		};

    		self.for_each_wrapped_interaction(func, get_indices, process_wrapped);
		}
    private:
    	LinkedCellsBatch<AsymBatch, SymBatch> batch;

    };
}


namespace april::container {
    struct LinkedCellsSoA : internal::LinkedCellsConfig{
        using ConfigT = LinkedCellsSoA;

        template <class U>
        using impl = internal::LinkedCellsSoAImpl<ConfigT, U>;
    };
}