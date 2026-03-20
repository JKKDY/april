#pragma once

#include <omp.h>
#include "april/exec/executors/executor_traits.hpp"

namespace april::exec {
    struct OmpExecutor {
        template<IsWorkAtom F>
        void execute(const size_t batch_count, F&& task) const {
            #pragma omp parallel for schedule(guided) num_threads(4)
            for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                task(i);
            }
        }

        [[nodiscard]] size_t num_threads() const noexcept {
            return static_cast<size_t>(omp_get_max_threads());
        }
    };
}