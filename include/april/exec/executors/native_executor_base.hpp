#pragma once
#include <atomic>
#include <thread>
#include <vector>
#include <algorithm>
#include <type_traits>

#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"


#if defined(_WIN32)
    #include <windows.h>
#elif defined(__linux__) || defined(__APPLE__)
    #include <pthread.h>
#endif

namespace april::exec::internal {

    inline void pin_thread_to_core(std::thread::native_handle_type handle, int core_id) {
        #if defined(_WIN32)
            // Windows uses a bitmask where bit 0 is core 0, bit 1 is core 1, etc.
            DWORD_PTR mask = (1ull << core_id);
            SetThreadAffinityMask((HANDLE)handle, mask);

        #elif defined(__linux__)
            // Linux uses a CPU set structure
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(core_id, &cpuset);
            pthread_setaffinity_np(handle, sizeof(cpu_set_t), &cpuset);

        #elif defined(__APPLE__)
            // macOS thread pinning is notoriously convoluted (mach thread binding).
            // Apple Silicon is heavily optimized to ignore standard pinning anyway.
            // For macOS, it is standard practice to do nothing and let macOS handle it.
            (void)handle;
            (void)core_id;
        #endif
    }

    inline void pin_current_thread(const int core_id) {
        #if defined(_WIN32)
            pin_thread_to_core(GetCurrentThread(), core_id);
        #elif defined(__linux__) || defined(__APPLE__)
            pin_thread_to_core(pthread_self(), core_id);
        #endif
    }
}

namespace april::exec::internal {
    // cross platform no-op
    // x86 / x86_64 (Intel & AMD)
    #if defined(__x86_64__) || defined(__i386__)
        #include <emmintrin.h>
        inline void cpu_pause() noexcept {
            _mm_pause(); // should keep power consumption low. see https://coffeebeforearch.github.io/2020/11/07/spinlocks-4.html
        }
        // ARM / AArch64 (Apple Silicon, AWS Graviton)
    #elif defined(__aarch64__) || defined(__arm__)
        inline void cpu_pause() noexcept {
            __asm__ volatile("yield" ::: "memory");
        }

        // Fallback for unknown architectures
    #else
        inline void cpu_pause() noexcept {
            // Compiler memory barrier to prevent the spin-loop
            // from being aggressively optimized away.
            __asm__ volatile("" ::: "memory");
        }
    #endif

} // namespace april::exec::internal

namespace april::exec::internal {

    struct TaskWrapper {
        void* ctx = nullptr;
        void (*invoke)(void*, size_t) = nullptr;
    };

    class NativeExecutorBase {
    public:
        [[nodiscard]] size_t num_threads() const noexcept { return threads.size() + 1; }

    protected:
        std::vector<std::jthread> threads;
        std::atomic<bool> terminate{false};

        // Mutable read-only state for current task
        mutable size_t total_tasks{0};
        mutable size_t chunk_size{1};
        mutable TaskWrapper active_task{};

        // Mutable contended state (Isolated in its own cache line)
        alignas(CACHE_LINE_SIZE) mutable std::atomic<size_t> current_idx{0};

        // Sets up the type-erased lambda and calculates chunk sizes
        template <IsWorkAtom F>
        void prepare_task(const size_t batch_count, F&& task) const {
            active_task.ctx = &task; // type erase the task
            active_task.invoke = [](void* ctx, size_t i) {
                // map type erased task back to original type and call with index i
                (*static_cast<std::remove_reference_t<F>*>(ctx))(i);
            };

            total_tasks = batch_count;
            current_idx.store(0, std::memory_order_relaxed);

            // TODO add more parameters to heuristic chunking (let users decide behaviour)
            chunk_size = std::max<size_t>(1, batch_count / ((threads.size() + 1) * 4));
        }

        // Shared work loop for both workers and the main thread.
        void process_tasks() const {
            while (true) {
                // fetch work until queue is empty
                const size_t start = current_idx.fetch_add(chunk_size, std::memory_order_relaxed);
                if (start >= total_tasks) break;

                // process fetched work
                const size_t end = std::min(start + chunk_size, total_tasks);
                for (size_t i = start; i < end; ++i) {
                    active_task.invoke(active_task.ctx, i);
                }
            }
        }
    };

} // namespace april::exec::internal