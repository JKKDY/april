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
		// The following hilbert curve encoding is based on Princeton numpy-hilbert-curve library:
		//	https://github.com/PrincetonLIPS/numpy-hilbert-curve
		// which in turn is based on Skillings grey-code "correction" procedure presented in:
		//	Skilling, J. (2004, April). Programming the Hilbert curve. In AIP Conference Proceedings
		//	(Vol. 707, No. 1, pp. 381-387). American Institute of Physics.

		// given a point and size of the grid (2^num_bits), return the hilbert index of that point
		// (num steps along the hilbert curve)
		inline uint64_t hilbert_encode(std::vector<uint32_t> coords, const int num_bits) {
		    const size_t num_dims = coords.size();

		    // Safety check
		    if (num_dims * num_bits > 64) {
		        throw std::overflow_error("Hilbert index would exceed 64 bits.");
		    }

			// Loop through all bits starting at the MSB (Most Significant Bit).
			// This is akin to iterating through the most coarse (hyper-) quadrant
			// and zooming in to divide it in each step.
			for (int i = num_bits - 1; i >= 0; --i) {

				// lower_mask = 0...001...1 with i least-significant bits set to 1.
				// Using this in an XOR operation (x ^ lower_mask) will flip all i
				// least significant bits (effectively processing the sub-quadrants).
				const uint32_t lower_mask = (1U << i) - 1;

				// Loop through all dimensions.
				// For each quadrant, we need to "unrotate" (or "unflip") its orientation
				// so the next iteration can treat it as a standard non-rotated box.
				for (size_t d = 0; d < num_dims; ++d) {

					// Check if the coordinate along the d-th dimension has its i-th bit set.
					// In other words: Is the point in the "Top" half of the current square?
					if ((coords[d] >> i) & 1) {

						// Yes: perform a horizontal reflection (invert the i-th lowest bits).
						// Intuition: This typically happens in "exit" regions. We flip the
						// coordinates so the exit point aligns with the entry of the next block.
						coords[0] ^= lower_mask;

					} else {

						// No: perform a geometric Transpose (Swap axes) to rotate the frame.
						// We swap the lower bits of the Primary Axis (0) and Current Axis (d).
						// This handles the "entry" regions where the curve winds inwards.
						const uint32_t t = (coords[0] ^ coords[d]) & lower_mask;
						coords[0] ^= t;
						coords[d] ^= t;
					}
				}
			}

			// build grey index (hilbert index of point but in grey number system) from modified coordinates
		    uint64_t gray_index = 0;
			// loop through all bits starting at the most significant (msb; most left; determines the top level quadrant)
		    for (int i = num_bits - 1; i >= 0; --i) {
		    	// loop through all dimensions
		        for (size_t d = 0; d < num_dims; ++d) {
		            const uint64_t bit = (coords[d] >> i) & 1;
		            gray_index = (gray_index << 1) | bit;
		        }
		    }

		    // convert grey index back to binary
			uint64_t bin = gray_index;
			bin ^= bin >> 1;
			bin ^= bin >> 2;
			bin ^= bin >> 4;
			bin ^= bin >> 8;
			bin ^= bin >> 16;
			bin ^= bin >> 32;
			return bin;
		}

		// wrapper for 3d hilbert encode
		inline uint64_t hilbert_encode_3d(uint32_t x, uint32_t y, uint32_t z, const int depth) {
		    return hilbert_encode({x, y, z}, depth);
		}
	} // namespace internal


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
					const uint64_t hilbert_curve_index = internal::hilbert_encode_3d(x, y, z, bits);
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
