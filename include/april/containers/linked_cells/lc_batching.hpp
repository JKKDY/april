#pragma once

#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/containers/batching/common.hpp"


namespace april::container::internal {

	template<typename AsymmetricBatch, typename SymmetricBatch>
	struct LinkedCellsBatch : batching::BatchBase<exec::internal::ParallelTrait::None,
		AsymmetricBatch::vector_trait & SymmetricBatch::vector_trait> {

		template<ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
		void for_each_pair (Func && f) const {
			for (const auto& chunk : sym_chunks)
				chunk.template for_each_pair<P, E>(f);

			for (const auto & chunk : asym_chunks)
				chunk.template for_each_pair<P, E>(f);
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













