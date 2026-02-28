#pragma once

#include "april/exec/policy.hpp"

#include "april/containers/batching/common.hpp"


namespace april::container::internal {

	template<typename AsymmetricBatch, typename SymmetricBatch>
	struct LinkedCellsBatch : batching::BatchBase<exec::internal::ParallelTrait::None,
		exec::internal::VectorTrait::ScalarOnly | exec::internal::VectorTrait::VectorOnly> {

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, typename Func>
		void for_each_pair (Func && f) const {
			for (const auto& chunk : sym_chunks)
				chunk.template for_each_pair<Mask, P, E>(f);

			for (const auto & chunk : asym_chunks)
				chunk.template for_each_pair<Mask, P, E>(f);
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













