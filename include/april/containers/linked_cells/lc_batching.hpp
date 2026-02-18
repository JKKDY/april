#pragma once

#include "april/base/policy.hpp"
#include "april/particle/fields.hpp"
#include "april/containers/batching/common.hpp"


namespace april::container::internal {

	template<typename AsymmetricBatch, typename SymmetricBatch>
	struct LinkedCellsBatch : SerialBatch {

		template<ParticleField Mask, ParallelPolicy P, VectorPolicy V, typename Func>
		void for_each_pair (Func && f) const {
			for (const auto& chunk : sym_chunks)
				chunk.template for_each_pair<Mask, P, V>(f);

			for (const auto & chunk : asym_chunks)
				chunk.template for_each_pair<Mask, P, V>(f);
		}

		void clear() {
			sym_chunks.clear();
			asym_chunks.clear();
		}

		[[nodiscard]] bool empty() const {
			return sym_chunks.empty() && asym_chunks.empty();
		}

		std::vector<SymmetricBatch> sym_chunks;
		std::vector<AsymmetricBatch> asym_chunks;
	};
}




