#pragma once

#include <omp.h>

#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"

namespace april::exec {
    struct OmpExecutor {

        struct Config {
            size_t n_threads = N_CPU_THREADS;
        };

        explicit OmpExecutor(const Config & config):
            n_threads(config.n_threads) {}

        template<ParallelPolicy P=ParallelPolicy::Threaded, IsWorkAtom F>
        void execute(const size_t batch_count, F&& task) const {
            if constexpr (P == ParallelPolicy::Serial) {
                for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                    task(i);
                }
            } else {
                #pragma omp parallel for schedule(guided) num_threads(n_threads)
                for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                    task(i);
                }
            }
        }

        [[nodiscard]] size_t num_threads() const noexcept {
            return n_threads;
        }

    private:
        unsigned n_threads;
    };
}