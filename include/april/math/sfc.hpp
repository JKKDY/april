#pragma once

#include <cstdint>
#include <vector>
#include <stdexcept>

#include "april/math/vec3.hpp"

namespace april::math::sfc {

    // -----------------------------
    // 1. MORTON (Z-Curve) UTILITIES
    // -----------------------------


    [[nodiscard]] AP_FORCE_INLINE uint64_t split_by_3(const uint32_t a) {
        uint64_t x = a & 0x1fffff; // mask to 21 bits
        x = (x | x << 32) & 0x1f00000000ffff;
        x = (x | x << 16) & 0x1f0000ff0000ff;
        x = (x | x <<  8) & 0x100f00f00f00f00f;
        x = (x | x <<  4) & 0x10c30c30c30c30c3;
        x = (x | x <<  2) & 0x1249249249249249;
        return x;
    }

    [[nodiscard]] AP_FORCE_INLINE uint64_t morton_key(const uint32_t x, const uint32_t y, const uint32_t z) {
        return split_by_3(x) | (split_by_3(y) << 1) | (split_by_3(z) << 2);
    }

    [[nodiscard]] AP_FORCE_INLINE uint64_t morton_key(const Vec3<uint32_t>& v) {
        return morton_key(v.x, v.y, v.z);
    }


    // --------------------------
    // 2. HILBERT CURVE UTILITIES
    // --------------------------
    // The following hilbert curve encoding is based on Princeton numpy-hilbert-curve library:
    // https://github.com/PrincetonLIPS/numpy-hilbert-curve
    // which in turn is based on Skillings grey-code "correction" procedure presented in:
    // Skilling, J. (2004, April). Programming the Hilbert curve. In AIP Conference Proceedings
    // (Vol. 707, No. 1, pp. 381-387). American Institute of Physics.

    // given a point and size of the grid (2^num_bits), return the hilbert index of that point
    // (num steps along the hilbert curve)
    // Optimized to take Vec3<T> to avoid std::vector heap allocations.
    template <typename T>
    [[nodiscard]] AP_FORCE_INLINE uint64_t hilbert_key(Vec3<T> coords, const int num_bits) {
        static_assert(std::is_integral_v<T>, "Hilbert coordinates must be integral");

        constexpr size_t num_dims = 3;

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
            const T lower_mask = (1U << i) - 1;

            // Loop through all dimensions.
            // For each quadrant, we need to "unrotate" (or "unflip") its orientation
            // so the next iteration can treat it as a standard non-rotated box.
            for (size_t d = 0; d < num_dims; ++d) {

                // Check if the coordinate along the d-th dimension has its i-th bit set.
                // In other words: Is the point in the "Top" half of the current square?
                if ((coords[static_cast<int>(d)] >> i) & 1) {

                    // Yes: perform a horizontal reflection (invert the i-th lowest bits).
                    // Intuition: This typically happens in "exit" regions. We flip the
                    // coordinates so the exit point aligns with the entry of the next block.
                    coords[0] ^= lower_mask;

                } else {

                    // No: perform a geometric Transpose (Swap axes) to rotate the frame.
                    // We swap the lower bits of the Primary Axis (0) and Current Axis (d).
                    // This handles the "entry" regions where the curve winds inwards.
                    const T t = (coords[0] ^ coords[static_cast<int>(d)]) & lower_mask;
                    coords[0] ^= t;
                    coords[static_cast<int>(d)] ^= t;
                }
            }
        }

        // build grey index (hilbert index of point but in grey number system) from modified coordinates
        uint64_t gray_index = 0;
        // loop through all bits starting at the most significant (msb; most left; determines the top level quadrant)
        for (int i = num_bits - 1; i >= 0; --i) {
            // loop through all dimensions
            for (size_t d = 0; d < num_dims; ++d) {
                const uint64_t bit = (coords[static_cast<int>(d)] >> i) & 1;
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
} // namespace april::math::sfc
