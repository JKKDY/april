#pragma once
#include <atomic>
#include <thread>
#include <vector>
#include <stdexcept>

#include "april/exec/threading/threading_context.hpp"
#include "april/exec/hardware.hpp"
#include "april/exec/threading/executor_concepts.hpp"
#include "internal/native_executor_base.hpp"
#include "april/exec/policy.hpp"


namespace april::exec {
    class NativeSpinExecutor : public internal::NativeExecutorBase {
    public:
        struct Config {
            size_t n_threads = default_thread_count;
            bool pin_threads = true;
        };

        explicit NativeSpinExecutor(const Config & config) {
            if (config.n_threads < 1) {
                throw std::invalid_argument("Executor thread count must be greater than zero.");
            }

            if (config.pin_threads) internal::pin_current_thread(0);
            threads.reserve(config.n_threads - 1);

            for (unsigned thread_idx = 1; thread_idx < config.n_threads; ++thread_idx) {
                threads.emplace_back(&NativeSpinExecutor::worker_loop, this, thread_idx);
                if (config.pin_threads) internal::pin_thread_to_core(threads.back().native_handle(), thread_idx);
            }
        }

        ~NativeSpinExecutor() {
            // set termination flag and drop the spin-lock barrier so sleeping threads wake up and see the flag
            terminate.store(true, std::memory_order_relaxed);
            run_signal.fetch_add(1, std::memory_order_release);
        }

        template <ParallelPolicy P = ParallelPolicy::Threaded, IsIndexedWork F>
        void execute(const size_t batch_count, F&& task) const {
            if (batch_count == 0) return;
            if constexpr (P == ParallelPolicy::Serial) {
                for (int i = 0; i < static_cast<int>(batch_count); ++i) {
                    task(i);
                }
            } else {

                prepare_task(batch_count, std::forward<F>(task));

                // reset completion counter, wake up threads and start processing
                threads_finished.store(0, std::memory_order_relaxed);
                run_signal.fetch_add(1, std::memory_order_release);

                internal::ScopedThreadContext ctx(0);
                process_tasks();

                // spin until all threads have finished
                const uint32_t target = static_cast<uint32_t>(threads.size());
                while (threads_finished.load(std::memory_order_acquire) != target) {
                    internal::cpu_pause();
                }
            }
        }

    private:
        // synchronization signals
        alignas(assumed_cache_line_size) mutable std::atomic<uint32_t> run_signal{0}; // increment to wake up threads
        alignas(assumed_cache_line_size) mutable std::atomic<uint32_t> threads_finished{0}; // if == #threads we are done

        void worker_loop(const int thread_idx) {
            internal::ScopedThreadContext ctx(thread_idx);
            uint32_t local_signal = 0; // used to compare to atomic (global) run signal. Only if they differ we run

            while (true) {
                // spin until signaled to start work again. If terminate signal received exit thread
                while (run_signal.load(std::memory_order_acquire) == local_signal) {
                    if (terminate.load(std::memory_order_relaxed)) return;
                    internal::cpu_pause();
                }
                local_signal = run_signal.load(std::memory_order_relaxed);

                // process tasks until queue empty. then signal completion.
                process_tasks();
                threads_finished.fetch_add(1, std::memory_order_release);
            }
        }
    };
} // namespace april::exec
