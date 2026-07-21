#pragma once

#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/kernel.hpp"
#include "april/containers/batching/batch.hpp"


namespace april::container::internal {

	template<typename AsymmetricBatch, typename SymmetricBatch>
	struct LinkedCellsBatch : batching::BatchBase<2,
		exec::execution_paths_intersection_t<
			typename AsymmetricBatch::execution_paths,
			typename SymmetricBatch::execution_paths
		>
	> {

		template<exec::ExecutionMode E, exec::IsKernel Func>
		void for_each (Func && f) const {
			for (const auto& chunk : sym_chunks)
				chunk.template for_each<E>(f);

			for (const auto & chunk : asym_chunks)
				chunk.template for_each<E>(f);
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















