#pragma once
#include <atomic>
#include <barrier>
#include <thread>
#include <vector>
#include <type_traits>
#include <algorithm>


#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"

namespace april::exec {
    class NativeBarrierExecutor {
    public:
        // n-1 workers because the main thread does work too
        explicit NativeBarrierExecutor(const unsigned int n = N_CPU_THREADS)
            : start_sync(n), end_sync(n) {
            threads.reserve(n-1);
            for (unsigned int i = 0; i < n-1; ++i) {
                threads.emplace_back(&NativeBarrierExecutor::worker_loop, this);
            }
        }

        ~NativeBarrierExecutor() {
            terminate.store(true, std::memory_order_relaxed);
            start_sync.arrive_and_wait();
        }


        template <IsWorkAtom F>
        void execute(const size_t batch_count, F&& task) const {
            if (batch_count == 0) return;

            // use type erasure to schedule task
            active_task.ctx = &task;
            active_task.invoke = [](void* ctx, size_t i) {
                (*static_cast<std::remove_reference_t<F>*>(ctx))(i);
            };

            total_tasks = batch_count;
            current_idx.store(0, std::memory_order_relaxed);

            // how many tasks each threads fetches at once
            chunk_size = std::max<size_t>(1, batch_count / ((threads.size() + 1) * 4));

            start_sync.arrive_and_wait(); // Wake workers
            process_tasks(); // Main thread processes tasks alongside workers
            end_sync.arrive_and_wait(); // Sync at the end of the phase
        }

        [[nodiscard]] size_t num_threads() const noexcept { return threads.size() + 1; }

    private:
        struct TaskWrapper {
            void* ctx = nullptr;
            void (*invoke)(void*, size_t) = nullptr;
        };

        std::vector<std::jthread> threads;
        std::atomic<bool> terminate{false};

        // mutable read-only state for current task
        mutable size_t total_tasks{0};
        mutable size_t chunk_size{1};
        mutable TaskWrapper active_task{};

        // mutable contended state (Isolated in its own cache line)
        alignas(CACHE_LINE_SIZE) mutable std::atomic<size_t> current_idx{0};

        // mutable barriers
        alignas(CACHE_LINE_SIZE) mutable std::barrier<> start_sync;
        alignas(CACHE_LINE_SIZE) mutable std::barrier<> end_sync;


        // The shared work loop for both workers and the main thread
        void process_tasks() const {
            while (true) {
                const size_t start = current_idx.fetch_add(chunk_size, std::memory_order_relaxed);
                if (start >= total_tasks) break;

                const size_t end = std::min(start + chunk_size, total_tasks);
                for (size_t i = start; i < end; ++i) {
                    active_task.invoke(active_task.ctx, i);
                }
            }
        }

        // thread function. processes tasks, then sleeps
        void worker_loop() const {
            while (true) {
                start_sync.arrive_and_wait();
                if (terminate.load(std::memory_order_relaxed)) return;

                process_tasks(); // Execute chunks

                end_sync.arrive_and_wait();
            }
        }
    };
} // namespace april::exec
