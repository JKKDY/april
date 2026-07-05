#pragma once
#include <atomic>
#include <thread>
#include <vector>
#include <algorithm>
#include <type_traits>

#include "april/exec/info.hpp"
#include "april/exec/executors/executor_traits.hpp"


#if defined(_WIN32)
#ifndef NOMINMAX
#  define NOMINMAX
#endif
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
        // Explicitly cast to native_handle_type to satisfy MinGW's pthread implementation
        pin_thread_to_core(
            reinterpret_cast<std::thread::native_handle_type>(GetCurrentThread()),
            core_id
        );
        #elif defined(__linux__) || defined(__APPLE__)
            pin_thread_to_core(pthread_self(), core_id);
        #endif
    }
}

namespace april::exec::internal {
    // cross platform no-op
    // x86 / x86_64 (Intel & AMD)
    #if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
        #if !defined(_MSC_VER)
            #include <immintrin.h>
        #endif
        inline void cpu_pause() noexcept {
            _mm_pause();
        }

    // ARM / AArch64 (Apple Silicon, AWS Graviton, Windows ARM)
    #elif defined(__aarch64__) || defined(__arm__) || defined(_M_ARM64) || defined(_M_ARM)
        inline void cpu_pause() noexcept {
            #if defined(_MSC_VER)
                __yield(); // MSVC ARM intrinsic
            #else
                __asm__ volatile("yield" ::: "memory"); // GCC/Clang ARM asm
            #endif
        }

    // Fallback for unknown architectures
    #else
        inline void cpu_pause() noexcept {
            #if defined(_MSC_VER)
                _ReadWriteBarrier(); // MSVC compiler memory barrier
            #else
                __asm__ volatile("" ::: "memory"); // GCC/Clang compiler memory barrier
            #endif
        }
    #endif
} // namespace april::exec::internal

namespace april::exec::internal {

    // apple clang has problems with jthread so well create our own
    struct worker_thread : std::thread {
        using std::thread::thread;
        ~worker_thread() { if (this->joinable()) this->join(); }
        worker_thread(worker_thread&&) noexcept = default;
        worker_thread& operator=(worker_thread&&) noexcept = default;
    };

    struct TaskWrapper {
        void* callable = nullptr; // need to store the task here because
        void (*invoke)(void*, size_t) = nullptr;
    };

    class NativeExecutorBase {
    public:
        [[nodiscard]] size_t num_threads() const noexcept { return threads.size() + 1; }

    protected:
        std::vector<worker_thread> threads;
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
            active_task.callable = &task; // type erase the task
            active_task.invoke = [](void* ctx, size_t task_idx) {
                // map type erased task back to original type and call with thread index and task index
                (*static_cast<std::remove_reference_t<F>*>(ctx))(task_idx);
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
                    active_task.invoke(active_task.callable, i);
                }
            }
        }
    };

} // namespace april::exec::internal