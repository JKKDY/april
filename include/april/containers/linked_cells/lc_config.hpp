#pragma once
#include <functional>
#include <optional>

#include "april/base/types.hpp"

namespace april::container {
	enum class CellSize {
		Cutoff,     // 1.0 * rc
		Half,       // 0.5 * rc
		Third,      // 0.33 * rc
		ManualAbs,  // Custom value (absolute)
		ManualFac   // Custom value (factor * rc)
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

		auto with_abs_cell_size(this auto&& self, const double cell_size) {
			self.manual_cell_size = cell_size;
			self.cell_size_strategy = CellSize::ManualAbs;
			return self;
		}

		auto with_cell_size_factor(this auto&& self, const double factor) {
			self.manual_cell_size = factor;
			self.cell_size_strategy = CellSize::ManualFac;
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
			case CellSize::ManualAbs: return manual_cell_size.value();
			case CellSize::ManualFac: return manual_cell_size.value() * rc;
			default: std::unreachable();
			}
		}
	};
}
