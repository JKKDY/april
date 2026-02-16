#pragma once
#include <cstdint>
#include "april/base/bitmask.hpp"

namespace april {
    enum class ComputePolicy : uint8_t {
         SerialScalar   = 0,
         ParallelScalar = 1 << 0,
         SerialVector   = 1 << 1,
         ParallelVector = (1 << 0) | (1 << 1), // Both bits set

         // Aliases for checking
         Parallel       = 1 << 0,
         Vector         = 1 << 1
    };

}




