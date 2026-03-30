/**
* @file scheduling.hpp
 * @brief Domain decomposition and thread-coloring schedulers.
 *
 * This file provides spatial coloring algorithms used to prevent data races
 * during parallel force accumulation with Newton's Third Law (N3L).
 *
 * By assigning cells to different "colors", we ensure that no two threads
 * simultaneously update the same particle or memory location.
 */

#pragma once
#include <stdexcept>
#include "april/base/types.hpp"


namespace april {

    /**
     * @brief 1-Color (Serial) scheduler.
     *
     * No parallel decomposition — all cells are processed in a single batch.
     * Zero synchronization overhead. Only safe when multi-threading is disabled.
     */
    inline size_t C01_schedule(const size_t /*x*/, const size_t /*y*/, const size_t /*z*/, const uint3& /*block_dim*/) {
        return 0;
    }

    /**
     * @brief 8-Color (Standard Parallel) scheduler.
     *
     * Classic 3D checkerboard coloring.
     *
     * @note Requires block_dim >= {2, 2, 2} to guarantee that the N3L stencil
     *       (radius 1 cell) does not cause overlapping updates on the same color.
     *
     * @throws std::invalid_argument if block dimensions are too small.
     */
    inline size_t C08_schedule(const size_t x, const size_t y, const size_t z, const uint3& block_dim) {
        if (block_dim.x < 2 || block_dim.y < 2 || block_dim.z < 2) {
            throw std::invalid_argument("[April] C08 scheduling requires block_size >= {2, 2, 2} to prevent data races.");
        }
        return (x % 2) + ((y % 2) << 1) + ((z % 2) << 2);
    }

    /**
     * @brief 18-Color (Dense Parallel) scheduler.
     *
     * Optimized for small domains or 1x1x1 blocks.
     * Uses a 3x3x2 modulo pattern to safely handle negative offsets in X and Y.
     */
    inline size_t C18_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 3) + (y % 3) * 3 + (z % 2) * 9;
    }

    /**
     * @brief 27-Color scheduler (for CellSize::Half).
     *
     * Designed for stencils that reach up to 2 cells in each direction.
     * Uses a 3x3x3 modulo pattern to guarantee race-free execution.
     */
    inline size_t C27_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 3) + (y % 3) * 3 + (z % 3) * 9;
    }

    /**
     * @brief 64-Color scheduler (for CellSize::Third).
     *
     * Intended for long-range stencils reaching up to 3 cells.
     * Uses a 4x4x4 modulo pattern for high-density neighbor searches.
     */
    inline size_t C64_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 4) + ((y % 4) << 2) + ((z % 4) << 4);
    }

    /**
     * @brief 2-Color "Slice" decomposition along the Z-axis.
     *
     * Decomposes the domain into XY-planes. Low synchronization overhead (only 2 barriers).
     * Well-suited for thin films or quasi-2D simulations.
     *
     * @throws std::invalid_argument if block_dim.z < 2.
     */
    inline size_t C02_Z_schedule(const size_t /*x*/, const size_t /*y*/, const size_t z, const uint3& block_dim) {
        if (block_dim.z < 2) {
            throw std::invalid_argument("[April] C02_Z scheduling requires block_size.z >= 2.");
        }
        return z % 2;
    }

    /**
     * @brief 4-Color "Pencil" decomposition along X and Y.
     *
     * Decomposes the domain into Z-columns (pencils).
     * Efficient for elongated geometries (e.g. nanowires, tubes) where Z >> X,Y.
     *
     * @throws std::invalid_argument if block_dim.x < 2 or block_dim.y < 2.
     */
    inline size_t C04_XY_schedule(const size_t x, const size_t y, const size_t /*z*/, const uint3& block_dim) {
        if (block_dim.x < 2 || block_dim.y < 2) {
            throw std::invalid_argument("[April] C04_XY scheduling requires block_size.x >= 2 and block_size.y >= 2.");
        }
        return (x % 2) + ((y % 2) << 1);
    }
}
