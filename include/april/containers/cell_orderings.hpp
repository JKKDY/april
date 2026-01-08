// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Julian Deller-Yee
#pragma once

#include <cstdint>

namespace april::container::internal {
	inline uint64_t split_by_3(const uint32_t a) {
		uint64_t x = a & 0x1fffff; // mask to 21 bits
		x = (x | x << 32) & 0x1f00000000ffff;
		x = (x | x << 16) & 0x1f0000ff0000ff;
		x = (x | x <<  8) & 0x100f00f00f00f00f;
		x = (x | x <<  4) & 0x10c30c30c30c30c3;
		x = (x | x <<  2) & 0x1249249249249249;
		return x;
	}

	inline uint64_t morton_3d_64(const uint32_t x, const uint32_t y, const uint32_t z) {
		return split_by_3(x) | (split_by_3(y) << 1) | (split_by_3(z) << 2);
	}

	void init_morton_order() {
		const size_t NX = cells_per_axis.x;
		const size_t NY = cells_per_axis.y;
		const size_t NZ = cells_per_axis.z;

		grid_to_linear_idx.resize(NX * NY * NZ);

		// 1. Create a temporary list of all coordinates
		struct Coord { size_t x, y, z; size_t original_flat_index; uint64_t morton; };
		std::vector<Coord> coords;
		coords.reserve(NX * NY * NZ);

		for (size_t z = 0; z < NZ; ++z) {
			for (size_t y = 0; y < NY; ++y) {
				for (size_t x = 0; x < NX; ++x) {
					// Calculate raw Morton code (the sparse one)
					uint64_t m = morton_3d_64(x, y, z);
					size_t flat = z * NY * NX + y * NX + x;
					coords.push_back({x, y, z, flat, m});
				}
			}
		}

		// 2. Sort based on Morton code
		// This removes the "gaps" because we are just sorting the existing valid cells
		std::sort(coords.begin(), coords.end(), [](const Coord& a, const Coord& b) {
			return a.morton < b.morton;
		});

		// 3. Build the lookup table
		// coords[i] is the i-th cell in memory (Sorted order).
		// We want to know: "Where in the sorted array is grid cell (x,y,z)?"
		for (size_t i = 0; i < coords.size(); ++i) {
			grid_to_linear_idx[coords[i].original_flat_index] = i;
		}
	}

	[[nodiscard]] uint32_t cell_pos_to_idx(this const auto & self, const uint32_t x, const uint32_t y, const uint32_t z) const noexcept{
		// get flat cell index
		uint32_t flat_idx = z * cells_per_axis.x * cells_per_axis.y + y * cells_per_axis.x + x;
		// map to cell ordering index (if possible)
		return self.cell_ordering.empty() ? flat_idx : self.cell_ordering[flat_idx];
	}
}
