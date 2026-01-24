#pragma once
#include <functional>
#include <optional>

#include "april/containers/cell_orderings.hpp"
#include "april/containers/batching.hpp"
#include "april/common.hpp"




namespace april::container {
	enum class CellSize {
		Cutoff,     // 1.0 * rc
		Half,       // 0.5 * rc
		Third,      // 0.33 * rc
		Fourth,     // 0.25 * rc
		Manual      // Custom value
	};
}

namespace april::container::internal {

	// ---------------------------------
	// LINKED CELLS CONFIGURATION STRUCT
	// ---------------------------------
	struct LinkedCellsConfig {
		CellSize cell_size_strategy = CellSize::Cutoff;
		std::optional<double> manual_cell_size;

		std::optional<std::function<std::vector<uint32_t>(uint3)>> cell_ordering_fn;
		uint3 block_size = {2,2,2};

		auto with_cell_size(this auto&& self, const double cell_size) {
			self.manual_cell_size = cell_size;
			self.cell_size_strategy = CellSize::Manual;
			return self;
		}

		auto with_cell_size(this auto&& self, const CellSize cell_size_strategy) {
			self.cell_size_strategy = cell_size_strategy;
			return self;
		}

		auto with_cell_ordering(this auto&& self, const std::function<std::vector<uint32_t>(uint3)> & ordering) {
			self.cell_ordering_fn = ordering;
			return self;
		}

		auto with_block_size(this auto&& self, const uint3 block_size) {
			self.block_size = block_size;
			return self;
		}

		auto with_block_size(this auto&& self, const uint3::type x, const uint3::type y, const uint3::type z) {
			self.block_size = uint3{x,y,z};
			return self;
		}

		auto with_block_size(this auto&& self, const uint3::type size) {
			self.block_size = uint3{size,size,size};
			return self;
		}

		[[nodiscard]] double get_width(const double rc) const {
			switch (cell_size_strategy) {
			case CellSize::Cutoff: return rc;
			case CellSize::Half:   return rc / 2.0;
			case CellSize::Third:  return rc / 3.0;
			case CellSize::Fourth:  return rc / 4.0;
			case CellSize::Manual: return manual_cell_size.value();
			default:               return rc;
			}
		}
	};


	// -----------------------------
	// LINKED CELLS INTERNAL STRUCTS
	// -----------------------------
	using cell_index_t = uint32_t;

	enum CellWrapFlag : uint8_t {
		NO_WRAP = 0,
		WRAP_X=1,
		WRAP_Y=2,
		WRAP_Z=4,
	};

	struct CellPair {
		const cell_index_t c1 = {};
		const cell_index_t c2 = {};
	};

	struct WrappedCellPair {
		const cell_index_t c1 = {};
		const cell_index_t c2 = {};
		const CellWrapFlag force_wrap = {};
		const vec3 shift = {};
	};



	// -----
	// Batch
	// -----
	template<typename Container>
	struct LinkedCellsBatch : batching::SerialBatch {
		explicit LinkedCellsBatch(Container & container) : container(container) {}

		template<env::FieldMask Mask, typename Func>
		void for_each_pair (Func && f) const {
			for (const auto& chunk : sym_chunks)
				chunk.template for_each_pair<Mask>(f);

			for (const auto & chunk : asym_chunks)
				chunk.template for_each_pair<Mask>(f);
		}

		void clear() {
			sym_chunks.clear();
			asym_chunks.clear();
		}

		[[nodiscard]] bool empty() const {
			return sym_chunks.empty() && asym_chunks.empty();
		}

		void add_sym_range(const batching::Range & range) {
			batching::SymmetricBatch<Container> batch (container);
			batch.range = range;
			sym_chunks.push_back(batch);
		}

		void add_asym_range(const batching::Range & range1, const batching::Range & range2) {
			batching::AsymmetricBatch<Container> batch (container);
			batch.range1 = range1;
			batch.range2 = range2;
			asym_chunks.push_back(batch);
		}

		std::vector<batching::SymmetricBatch<Container>> sym_chunks;
		std::vector<batching::AsymmetricBatch<Container>> asym_chunks;
	private:
		Container & container;
	};
}
