#pragma once

#include <cstdint>

namespace april {
    struct Layout {
        struct AoS {};
        struct SoA {};
        template<uint8_t ChunkSize = 8> struct AoSoA {
            static constexpr uint8_t chunk_size = ChunkSize;
        };
    };
}
