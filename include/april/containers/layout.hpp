#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>

#include "april/base/types.hpp"

namespace april {
    namespace internal {
        consteval size_t default_chunk_size() {
            // # elements in a simd vector
            constexpr size_t simd_width = packed::size();

            // # elements to fill a standard cache line
            constexpr size_t cache_line_elements = 64 / sizeof(packed::value_type);
            // static_assert(sizeof(packed::scalar_type) == 4);

            // chunk size must be at least a cache line
            return std::max(simd_width, cache_line_elements);
        }
    }

    /**
     * @brief Storage layout tags used to select a container memory layout.
     *
     * Layout tags are passed as template arguments to containers, for example:
     *
     * @code
     * auto container = LinkedCells<Layout::SoA>();
     * auto container = DirectSum<Layout::AoSoA<16>>();
     * @endcode
     */
    struct Layout {
        struct AoS {};
        struct SoA {};
        template<uint8_t ChunkSize = internal::default_chunk_size()> struct AoSoA {
            static constexpr uint8_t chunk_size = ChunkSize;
        };
    };
}


