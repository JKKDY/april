#pragma once
#include <stdexcept>
#include "april/base/types.hpp"


namespace april {

    // 1-Color (Serial)
    // No parallel scheduling. Only use if multi-threading is disabled.
    inline size_t C01_schedule(const size_t /*x*/, const size_t /*y*/, const size_t /*z*/, const uint3& /*block_dim*/) {
        return 0;
    }

    // 8-Color (Standard Parallel)
    // Requires a minimum block size of 2x2x2 to safely apply Newton's 3rd Law without atomic locks.
    inline size_t C08_schedule(const size_t x, const size_t y, const size_t z, const uint3& block_dim) {
        if (block_dim.x < 2 || block_dim.y < 2 || block_dim.z < 2) {
            throw std::invalid_argument("[April] C08 scheduling requires block_size >= {2, 2, 2} to prevent data races.");
        }
        return (x % 2) + ((y % 2) << 1) + ((z % 2) << 2);
    }

    // 18-Color (Dense Parallel)
    // Safe for 1x1x1 blocks when applying N3L. Accounts for negative offsets in X and Y.
    inline size_t C18_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 3) + (y % 3) * 3 + (z % 2) * 9;
    }

    // 27-Color (For CellSize::Half)
    // Stencil reaches 2 cells in each direction. Modulo 3 guarantees safety.
    inline size_t C27_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 3) + (y % 3) * 3 + (z % 3) * 9;
    }

    // 64-Color (For CellSize::Third)
    // Stencil reaches 3 cells in each direction. Modulo 4 guarantees safety.
    inline size_t C64_schedule(const size_t x, const size_t y, const size_t z, const uint3& /*block_dim*/) {
        return (x % 4) + ((y % 4) << 2) + ((z % 4) << 4);
    }

    // 2-Color "Slice" Decomposition (Z-axis only)
    // Great for thin films. Only 2 barriers per step. Threads lock whole XY planes.
    inline size_t C02_Z_schedule(const size_t /*x*/, const size_t /*y*/, const size_t z, const uint3& block_dim) {
        if (block_dim.z < 2) {
            throw std::invalid_argument("[April] C02_Z scheduling requires block_size.z >= 2.");
        }
        return z % 2;
    }

    // 4-Color "Pencil" Decomposition (X-Y only)
    // Great for long tubes/nanowires. 4 barriers per step. Threads process entire Z-columns.
    inline size_t C04_XY_schedule(const size_t x, const size_t y, const size_t /*z*/, const uint3& block_dim) {
        if (block_dim.x < 2 || block_dim.y < 2) {
            throw std::invalid_argument("[April] C04_XY scheduling requires block_size.x >= 2 and block_size.y >= 2.");
        }
        return (x % 2) + ((y % 2) << 1);
    }
}
