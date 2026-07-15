#pragma once

#include <omp.h>

#include "april/exec/hardware.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/threading/executor_concepts.hpp"
#include "april/exec/threading/threading_context.hpp"

namespace april::exec {
    struct OmpExecutor {

        struct Config {
            size_t n_threads = default_thread_count;
        };

        explicit OmpExecutor(const Config & config):
            n_threads(config.n_threads) {
            if (n_threads < 1) {
                throw std::invalid_argument("Executor thread count must be greater than zero.");
            }
        }

        template<ParallelPolicy P=ParallelPolicy::Threaded, IsIndexedWork F>
        void execute(const size_t batch_count, F&& task) const {
            if constexpr (P == ParallelPolicy::Serial) {
                for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                    task(i);
                }
            } else {
                #pragma omp parallel num_threads(n_threads)
                {
                    internal::ScopedThreadContext ctx(omp_get_thread_num());

                    #pragma omp for schedule(guided)
                    for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                        task(i);
                    }
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