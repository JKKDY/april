// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Julian Deller-Yee
#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "april/common.hpp"

namespace april::container {
	namespace internal {
		// --------------------------
		// MORTON (Z-Curve) UTILITIES
		// --------------------------
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


		// -----------------------
		// HILBERT CURVE UTILITIES
		// -----------------------
		// we use the following hilbert curve basis shape (fills a 2x2x2 octant)
		// Index (Time)   Coord (x,y,z)    Binary (val)      Geometric Octant Value (4z + 2y + x)
		//       0          (0, 0, 0)             000                           0
		//       1          (1, 0, 0)             001                           1
		//       2          (1, 1, 0)             011                           3
		//       3          (0, 1, 0)             010                           2
		//       4          (0, 1, 1)             110                           6
		//       5          (1, 1, 1)             111                           7
		//       6          (1, 0, 1)             101                           5
		//       7          (0, 0, 1)             100                           4
		// the following table the basis shape for every one of the possible 24 rotations
		// given a rotation, it gives the sequence of octants. These define a base hilbert curve.
		// the sequence is specified by the "time" (number of steps when traversing the curve) at which ocant
		// i (where i the geometric octant value as per the above grey code) should be visited
		// i.e. maps space -> time
		static constexpr uint8_t hilbert_octant_to_index[24][8] = {
			{0, 1, 3, 2, 7, 6, 4, 5}, {0, 7, 1, 6, 3, 4, 2, 5}, {0, 3, 7, 4, 1, 2, 6, 5},
			{2, 3, 1, 0, 5, 4, 6, 7}, {4, 3, 5, 2, 7, 0, 6, 1}, {6, 5, 1, 2, 7, 4, 0, 3},
			{4, 7, 5, 6, 3, 0, 2, 1}, {6, 7, 1, 0, 5, 4, 2, 3}, {0, 1, 3, 2, 7, 6, 4, 5},
			{2, 1, 3, 0, 5, 6, 4, 7}, {4, 5, 7, 6, 3, 2, 0, 1}, {6, 1, 7, 0, 5, 2, 4, 3},
			{0, 7, 1, 6, 3, 4, 2, 5}, {2, 1, 3, 0, 5, 6, 4, 7}, {4, 3, 5, 2, 7, 0, 6, 1},
			{4, 7, 5, 6, 3, 0, 2, 1}, {0, 1, 3, 2, 7, 6, 4, 5}, {0, 3, 7, 4, 1, 2, 6, 5},
			{2, 3, 1, 0, 5, 4, 6, 7}, {2, 1, 3, 0, 5, 6, 4, 7}, {4, 5, 7, 6, 3, 2, 0, 1},
			{4, 7, 5, 6, 3, 0, 2, 1}, {6, 5, 1, 2, 7, 4, 0, 3}, {6, 7, 1, 0, 5, 4, 2, 3}
		};

		// There are 24 orientation-preserving 90° rotations in 3D (octahedral group):
		// choose where the +X axis maps (6 choices), then choose a perpendicular +Y axis (4 choices) -> 6 x 4 = 24

		// consider the sub octants (cells) of a 2x2x2 block (octet). Each sub octant contains its own hilbert curve shape.
		// the octet contains its own hilbert curve that specifies which sub octant is adjacent to which.
		// sub octants are rotated such that the end of the previous matches the start of the next

		// the following is a state transition table for the 3D Hilbert curve (24 rotations x 8 octants)
		// given an orientation and a current sub octant, it gives the orientation of the next octant
		static constexpr uint8_t hilbert_rotation[24][8] = {
			{0,3,4,0,4,7,8,7}, {1,2,5,1,5,6,9,6}, {0,3,4,0,4,7,8,7}, {1,2,5,1,5,6,9,6},
			{10,13,14,10,14,17,18,17}, {11,12,15,11,15,16,19,16}, {10,13,14,10,14,17,18,17}, {11,12,15,11,15,16,19,16},
			{20,23,2,20,2,5,6,5}, {21,22,3,21,3,4,7,4}, {20,23,2,20,2,5,6,5}, {21,22,3,21,3,4,7,4},
			{0,3,4,0,4,7,8,7}, {1,2,5,1,5,6,9,6}, {0,3,4,0,4,7,8,7}, {1,2,5,1,5,6,9,6},
			{10,13,14,10,14,17,18,17}, {11,12,15,11,15,16,19,16}, {10,13,14,10,14,17,18,17}, {11,12,15,11,15,16,19,16},
			{20,23,2,20,2,5,6,5}, {21,22,3,21,3,4,7,4}, {20,23,2,20,2,5,6,5}, {21,22,3,21,3,4,7,4}
		};

		// calculates the Hilbert index for a point in a virtual 2^bits x 2^bits x 2^bits box.
		inline uint64_t hilbert_3d_64(const uint32_t x, const uint32_t y, const uint32_t z, const int bits) {
			uint64_t index = 0;
			int rotation = 0; // Canonical initial orientation

			// we build the index bit by bit
			// hilbert curves are self similar; we start at the most coarse curve and iteratively refine it (
			for (int i = bits - 1; i >= 0; --i) {
				const uint32_t mask = 1 << i; // get the most significant bit (most left)

				// which of the 8 octants are we in at the current (refinement) level?
				// first get the octant index (without considering rotations due to the curve)
				// if x has the bit set, it's on the "Right".
				// if y has the bit set, it's on the "Top".
				// if z has the bit set, it's in the "Front".
				const uint8_t octant = ((x & mask) ? 1 : 0) | ((y & mask) ? 2 : 0) | ((z & mask) ? 4 : 0);

				// get the octant index considering rotations
				const uint8_t seg_index = hilbert_octant_to_index[rotation][octant];

				// 3. Accumulate index
				index = (index << 3) | seg_index;

				// 4. Update rotation for the next level
				rotation = hilbert_rotation[rotation][octant];
			}
			return index;
		}
	}

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
					const uint64_t m = internal::morton_3d_64(x, y, z);
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
			cell_ordering[coords[i].original_flat_index] = i;
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
					const uint64_t hilbert_curve_index = internal::hilbert_3d_64(x, y, z, bits);
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
