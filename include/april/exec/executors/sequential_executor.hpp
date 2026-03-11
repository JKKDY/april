#pragma once

#include "april/exec/executors/executor_traits.hpp"

namespace april::exec {
    struct SequentialExecutor {
        template<IsWorkAtom F>
        void execute(const size_t batch_count, F&& task) const {
            for (size_t i = 0; i < batch_count; ++i) {
                task(i);
            }
        }

        [[nodiscard]] constexpr size_t num_threads() const noexcept {
            return 1;
        }
    };

} // namespace april::exec