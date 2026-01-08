// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Julian Deller-Yee
#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "april/common.hpp"

namespace april::container {
	namespace internal {
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
	}

	inline std::vector<uint32_t> morton_order(const uint3 & cells_per_axis) {
		const size_t N = cells_per_axis.x * cells_per_axis.y * cells_per_axis.z;
		std::vector<uint32_t> cell_ordering(N);

		// helper struct for storing, cell position, flat index, morton code
		struct Coord {
			uint32_t x, y, z;
			uint32_t original_flat_index;
			uint64_t morton;
		};

		// calculate flat index & morton code for every cell
		std::vector<Coord> coords;
		coords.reserve(N);
		for (uint32_t z = 0; z < cells_per_axis.z; ++z) {
			for (uint32_t y = 0; y <  cells_per_axis.y; ++y) {
				for (uint32_t x = 0; x <  cells_per_axis.x; ++x) {
					// Calculate raw Morton code (the sparse one)
					const uint64_t m = internal::morton_3d_64(x, y, z);
					const uint32_t flat = z * cells_per_axis.y * cells_per_axis.x + y * cells_per_axis.x + x;
					coords.push_back({x, y, z, flat, m});
				}
			}
		}

		// sort based on Morton code
		std::ranges::sort(coords, [](const Coord& a, const Coord& b) {
			return a.morton < b.morton;
		});

		// build the lookup table
		// coords[i] is the i-th cell in memory (Sorted order).
		// We want to know: "Where in the sorted array is grid cell (x,y,z)?"
		for (size_t i = 0; i < coords.size(); ++i) {
			cell_ordering[coords[i].original_flat_index] = i;
		}

		return cell_ordering;
	}

}
