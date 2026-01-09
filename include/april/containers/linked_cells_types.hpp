#pragma once
#include <functional>
#include <optional>

#include "april/containers/cell_orderings.hpp"
#include "april/containers/batching.hpp"
#include "april/common.hpp"

namespace april::container::internal {
	struct LinkedCellsConfig {
		std::optional<double> cell_size_hint;
		std::optional<std::function<std::vector<uint32_t>(uint3)>> cell_ordering_fn;
		uint8_t super_batch_size = 1;

		auto with_cell_size(this auto&& self, const double cell_size) {
			self.cell_size_hint = cell_size;
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
	};

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

	struct AsymmetricChunk {
		std::ranges::iota_view<size_t, size_t> indices1;
		std::ranges::iota_view<size_t, size_t> indices2;
	};

	struct SymmetricChunk {
		std::ranges::iota_view<size_t, size_t> indices;
	};

	struct AsymmetricChunkedBatch : SerialBatch<BatchSymmetry::Asymmetric> {
		std::vector<AsymmetricChunk> chunks;
	};

	struct SymmetricChunkedBatch : SerialBatch<BatchSymmetry::Symmetric> {
		std::vector<SymmetricChunk> chunks;
	};

	struct AsymmetricBatch : SerialBatch<BatchSymmetry::Asymmetric> {
		std::ranges::iota_view<size_t, size_t> indices1;
		std::ranges::iota_view<size_t, size_t> indices2;
	};

}
