#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>
#include <bit>

#include "april/base/types.hpp"
#include "april/math/sfc.hpp"

namespace april::container {
	// ----------------
	// MORTON (Z-Curve)
	// ----------------
	inline std::vector<uint32_t> morton_order(const uint3 & cells_per_axis) {
		const size_t N = cells_per_axis.x * cells_per_axis.y * cells_per_axis.z;
		std::vector<uint32_t> cell_ordering(N);

		// helper struct for storing flat index & morton code
		struct Coord {
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
					const uint64_t m = math::sfc::morton_key(x, y, z);
					const uint32_t flat = z * cells_per_axis.y * cells_per_axis.x + y * cells_per_axis.x + x;
					coords.push_back({flat, m});
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
			cell_ordering[coords[i].original_flat_index] = static_cast<uint32_t>(i);
		}

		return cell_ordering;
	}


	// -------------
	// HILBERT CURVE
	// -------------
	inline std::vector<uint32_t> hilbert_order(const uint3 & cells_per_axis) {
		const size_t N = cells_per_axis.x * cells_per_axis.y * cells_per_axis.z;
		std::vector<uint32_t> cell_ordering(N);

		// Determine necessary depth (Order of curve)
		const uint32_t max_dim = std::max({cells_per_axis.x, cells_per_axis.y, cells_per_axis.z});
		// bits = ceil(log2(max_dim)), but handling the 0-based index range
		const int bits = std::bit_width(max_dim > 0 ? max_dim - 1 : 0);

		// helper struct for storing flat index & hilbert curve index
		struct Coord {
			uint32_t original_flat_index;
			uint64_t hilbert_curve_index;
		};

		std::vector<Coord> coords;
		coords.reserve(N);

		for (uint32_t z = 0; z < cells_per_axis.z; ++z) {
			for (uint32_t y = 0; y < cells_per_axis.y; ++y) {
				for (uint32_t x = 0; x < cells_per_axis.x; ++x) {
					// Calculate Hilbert Key assuming the point is inside the virtual power-of-two box
					const uint64_t hilbert_curve_index = math::sfc::hilbert_key(uint3{x, y, z}, bits);
					const uint32_t flat = z * cells_per_axis.y * cells_per_axis.x + y * cells_per_axis.x + x;
					coords.push_back({flat, hilbert_curve_index});
				}
			}
		}

		// sort based on hilbert curve index
		std::ranges::sort(coords, [](const Coord& a, const Coord& b) {
		   return a.hilbert_curve_index < b.hilbert_curve_index;
		});

		// build the lookup table
		// coords[i] is the i-th cell in memory (Sorted order).
		// We want to know: "Where in the sorted array is grid cell (x,y,z)?"
		for (size_t i = 0; i < coords.size(); ++i) {
			cell_ordering[coords[i].original_flat_index] = static_cast<uint32_t>(i);
		}

		return cell_ordering;
	}

}
