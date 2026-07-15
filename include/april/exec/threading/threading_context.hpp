#pragma once


namespace april::exec {
    namespace internal {
        // The hidden thread local storage state for thread id
        inline thread_local int current_thread_id = -1;

        // RAII Guard to ensure safe setup and teardown
        struct ScopedThreadContext {
            int previous_id; // store previous id so we can restore original state. Needed in case of nesting
            explicit ScopedThreadContext(const int new_id) noexcept
                : previous_id(current_thread_id) {
                current_thread_id = new_id;
            }
            ~ScopedThreadContext() {
                current_thread_id = previous_id; // restore previous state
            }
            ScopedThreadContext(const ScopedThreadContext&) = delete;
            ScopedThreadContext& operator=(const ScopedThreadContext&) = delete;
        };
    }

    // get current index of thread (within its executor) (unsafe: may return -1. Mainly used for debugging)
    [[nodiscard]] inline int thread_index_direct() noexcept {
        return internal::current_thread_id;
    }

    // get current index of thread (within its executor) (safe)
    [[nodiscard]] inline unsigned thread_index() noexcept {
        return internal::current_thread_id < 0 ? 0 : internal::current_thread_id;
    }

    // check if thread is running parallel (within the executor)
    [[nodiscard]] inline bool is_parallel() noexcept {
        return internal::current_thread_id != -1;
    }
}