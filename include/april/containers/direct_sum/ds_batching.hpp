#pragma once

#include <vector>

namespace april::container::internal {

    template <typename SymmetricBatch, typename AsymmetricBatch>
    struct SymmetricTaskGroup {
        std::vector<SymmetricBatch> diagonals; // Phase 0
        std::vector<std::vector<AsymmetricBatch>> off_diagonals;  // Phases 1 to N-1
    };

    template <typename AsymmetricBatch>
    struct AsymmetricTaskGroup {
        std::vector<std::vector<AsymmetricBatch>> phases;
    };

}
