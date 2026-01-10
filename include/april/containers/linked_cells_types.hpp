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
		uint8_t super_batch_size = 1;

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

		auto with_super_batch_size(this auto&& self, const uint8_t super_batch_size) {
			self.super_batch_size = super_batch_size;
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


	// ------------------------------
	// Batch Work Units (The "Atoms")
	// ------------------------------
	struct AsymmetricChunk {
		std::ranges::iota_view<size_t, size_t> indices1;
		std::ranges::iota_view<size_t, size_t> indices2;
	};

	struct SymmetricChunk {
		std::ranges::iota_view<size_t, size_t> indices;
	};


	// --------------
	// Compound Batch
	// --------------
	struct UnifiedBatch : SerialBatch<BatchType::Compound> {
		std::pair<env::ParticleType, env::ParticleType> types;

		std::vector<AsymmetricChunk> sym_chunks;
		std::vector<SymmetricChunk> asym_chunks;

		void clear() {
			sym_chunks.clear();
			asym_chunks.clear();
		}

		[[nodiscard]] bool empty() const {
			return sym_chunks.empty() && asym_chunks.empty();
		}
	};

}
