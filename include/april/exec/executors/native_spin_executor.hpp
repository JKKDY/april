#pragma once
#include <atomic>
#include <thread>
#include <vector>

#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"
#include "native_executor_base.hpp"
#include "april/exec/policy.hpp"


namespace april::exec {
    class NativeSpinExecutor : public internal::NativeExecutorBase {
    public:
        struct Config {
            size_t n_threads = N_CPU_THREADS;
            bool pin_threads = true;
        };

        explicit NativeSpinExecutor(const Config & config) {
            if (config.pin_threads) internal::pin_current_thread(0);
            threads.reserve(config.n_threads - 1);

            for (unsigned int i = 0; i < config.n_threads - 1; ++i) {
                threads.emplace_back(&NativeSpinExecutor::worker_loop, this);
                if (config.pin_threads) internal::pin_thread_to_core(threads.back().native_handle(), i + 1);
            }
        }

        ~NativeSpinExecutor() {
            // set termination flag and drop the spin-lock barrier so sleeping threads wake up and see the flag
            terminate.store(true, std::memory_order_relaxed);
            run_signal.fetch_add(1, std::memory_order_release);
        }

        template <ParallelPolicy P = ParallelPolicy::Threaded, IsWorkAtom F>
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
        alignas(CACHE_LINE_SIZE) mutable std::atomic<uint32_t> run_signal{0}; // increment to wake up threads
        alignas(CACHE_LINE_SIZE) mutable std::atomic<uint32_t> threads_finished{0}; // if == #threads we are done

        void worker_loop() {
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
