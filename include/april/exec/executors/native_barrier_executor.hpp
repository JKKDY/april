#pragma once
#include <atomic>
#include <barrier>
#include <thread>
#include <vector>

#include "context.hpp"
#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"
#include "native_executor_base.hpp"
#include "april/exec/policy.hpp"


namespace april::exec {
    class NativeBarrierExecutor : public internal::NativeExecutorBase {
    public:
        struct Config {
            size_t n_threads = N_CPU_THREADS;
            bool pin_threads = true;
        };

        explicit NativeBarrierExecutor(const Config & config)
            : start_sync(config.n_threads), end_sync(config.n_threads) {

            if (config.pin_threads) internal::pin_current_thread(0);
            threads.reserve(config.n_threads - 1);

            for (unsigned int i = 0; i < config.n_threads - 1; ++i) {
                threads.emplace_back(&NativeBarrierExecutor::worker_loop, this, i+1);   // +1 because main thread is 0
                if (config.pin_threads) internal::pin_thread_to_core(threads.back().native_handle(), i + 1);
            }
        }

        ~NativeBarrierExecutor() {
            // signal termination
            terminate.store(true, std::memory_order_relaxed);
            start_sync.arrive_and_wait();
        }

        template <ParallelPolicy P=ParallelPolicy::Threaded, IsWorkAtom F>
        void execute(const size_t batch_count, F&& task) const {
            if (batch_count == 0) return;
            if constexpr (P == ParallelPolicy::Serial) {
                for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                    task(i);
                }
            } else {
                internal::ScopedThreadContext ctx(0);
                prepare_task(batch_count, std::forward<F>(task));
                start_sync.arrive_and_wait(); // Wake workers
                process_tasks();    // Main thread works as well
                end_sync.arrive_and_wait();   // Sync at phase end
            }
        }

    private:
        alignas(CACHE_LINE_SIZE) mutable std::barrier<> start_sync;
        alignas(CACHE_LINE_SIZE) mutable std::barrier<> end_sync;

        void worker_loop(const int thread_idx) const {
            internal::ScopedThreadContext ctx(thread_idx);
            while (true) {
                start_sync.arrive_and_wait();
                if (terminate.load(std::memory_order_relaxed)) return;

                process_tasks();

                end_sync.arrive_and_wait();
            }
        }
    };
} // namespace april::exec

