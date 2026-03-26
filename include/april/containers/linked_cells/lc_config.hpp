#pragma once
#include <functional>
#include <optional>

#include "april/base/types.hpp"
#include "lc_scheduling.hpp"

namespace april::container {
	enum class CellSize {
		Cutoff,     // 1.0 * rc
		Half,       // 0.5 * rc
		Third,      // 0.33 * rc
		ManualAbs,  // Custom value (absolute)
		ManualFac   // Custom value (factor * rc)
	};

	enum class SkinSize {
		Factor,
		Absolute
	 };
}

namespace april::container::internal {

	// ---------------------------------
	// LINKED CELLS CONFIGURATION STRUCT
	// ---------------------------------
	struct LinkedCellsConfig {
		CellSize cell_size_strategy = CellSize::Cutoff;
		std::optional<double> manual_cell_size;

		SkinSize skin_strategy = SkinSize::Factor;
		double skin_value = 0.1;

		std::optional<std::function<std::vector<uint32_t>(uint3)>> cell_ordering_fn;
		uint3 block_size = {2,2,2};

		std::function<size_t(size_t, size_t, size_t, uint3)> schedule_phases = C08_schedule;

		auto&& with_abs_cell_size(this auto&& self, const double cell_size) {
			self.manual_cell_size = cell_size;
			self.cell_size_strategy = CellSize::ManualAbs;
			return self;
		}

		auto&& with_cell_size_factor(this auto&& self, const double factor) {
			self.manual_cell_size = factor;
			self.cell_size_strategy = CellSize::ManualFac;
			return self;
		}

		auto&& with_cell_size(this auto&& self, const CellSize cell_size_strategy) {
			self.cell_size_strategy = cell_size_strategy;
			return self;
		}

		auto&& with_cell_ordering(this auto&& self, const std::function<std::vector<uint32_t>(uint3)> & ordering) {
			self.cell_ordering_fn = ordering;
			return self;
		}

		auto&& with_block_size(this auto&& self, const uint3 block_size) {
			self.block_size = block_size;
			return self;
		}

		auto&& with_block_size(this auto&& self, const uint3::type x, const uint3::type y, const uint3::type z) {
			self.block_size = uint3{x,y,z};
			return self;
		}

		auto&& with_block_size(this auto&& self, const uint3::type size) {
			self.block_size = uint3{size,size,size};
			return self;
		}

		auto&& with_skin_factor(this auto&& self, const double factor) {
			// factor is a multiplier of rc
			self.skin_value = factor;
			self.skin_strategy = SkinSize::Factor;
			return self;
		}

		auto&& with_absolute_skin(this auto&& self, const double skin) {
			// absolute distance buffer
			self.skin_value = skin;
			self.skin_strategy = SkinSize::Absolute;
			return self;
		}

		auto&& with_scheduling(this auto&& self, std::function<size_t(size_t, size_t, size_t)> fn) {
			self.schedule_phases = std::move(fn);
			return self;
		}

		[[nodiscard]] double get_width(const double max_force_cutoff) const {
			switch (cell_size_strategy) {
			case CellSize::Cutoff: return max_force_cutoff;
			case CellSize::Half:   return max_force_cutoff / 2.0;
			case CellSize::Third:  return max_force_cutoff / 3.0;
			case CellSize::ManualAbs: return manual_cell_size.value();
			case CellSize::ManualFac: return manual_cell_size.value() * max_force_cutoff;
			default: std::unreachable();
			}
		}

		[[nodiscard]] double get_skin(const double cell_size) const {
			switch (skin_strategy) {
			case SkinSize::Factor: return skin_value * cell_size;
			case SkinSize::Absolute: return skin_value;
			default: std::unreachable();
			}
		}
	};
}

