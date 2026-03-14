#pragma once

#include <vector>

namespace april::container::internal {

    template <typename SymmetricBatch, typename AsymmetricBatch>
    struct SymTaskGroup {
        std::vector<SymmetricBatch> diagonals; // Phase 0
        std::vector<std::vector<AsymmetricBatch>> off_diagonals;  // Phases 1 to N-1
    };

    template <typename AsymmetricBatch>
    struct AsymTaskGroup {
        std::vector<std::vector<AsymmetricBatch>> phases;
    };

}
