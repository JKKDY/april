#pragma once
#include <atomic>
#include <barrier>
#include <thread>
#include <vector>
#include <type_traits>
#include <algorithm>

namespace april::exec {
    class NativeExecutor {
    public:
        // n-1 workers because the main thread does work too
        explicit NativeExecutor(const unsigned int n = std::max<unsigned int>(1, n_threads - 1))
            : start_sync(n + 1), end_sync(n + 1) {
            threads.reserve(n);
            for (unsigned int i = 0; i < n; ++i) {
                threads.emplace_back(&NativeExecutor::worker_loop, this);
            }
        }

        ~NativeExecutor() {
            terminate.store(true, std::memory_order_relaxed);
            start_sync.arrive_and_wait();
        }

        template <typename F>
        void execute(const size_t batch_count, F&& task) {
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

        // read only state
        size_t total_tasks{0};
        size_t chunk_size{1};
        TaskWrapper active_task{};

        // contended state (Isolated in its own cache line)
        alignas(cache_line_size) std::atomic<size_t> current_idx{0};

        // barriers
        alignas(cache_line_size) std::barrier<> start_sync;
        alignas(cache_line_size) std::barrier<> end_sync;

        // The shared work loop for both workers and the main thread
        void process_tasks() {
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
        void worker_loop() {
            while (true) {
                start_sync.arrive_and_wait();
                if (terminate.load(std::memory_order_relaxed)) return;

                process_tasks(); // Execute chunks

                end_sync.arrive_and_wait();
            }
        }
    };
} // namespace april::exec
