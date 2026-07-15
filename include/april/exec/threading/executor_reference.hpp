#pragma once

#include <stddef.h>
#include "april/exec/policy.hpp"


namespace april::exec {

    template <typename RealExecutor>
    class ThreadExecutorRef {
    public:
        using Config =  RealExecutor::Config;

        ThreadExecutorRef() = default;
        explicit ThreadExecutorRef(RealExecutor* p) : ptr(p) {}

        void bind(RealExecutor* p) {
            ptr = p;
        }

        [[nodiscard]] size_t num_threads() const noexcept {
            APRIL_ASSERT(ptr != nullptr, "[APRIL] ExecutorRef executed before being bound to System!");
            return ptr->num_threads();
        }

        [[nodiscard]] constexpr bool is_bound() const noexcept {
            return ptr != nullptr;
        }


        template <ParallelPolicy P = ParallelPolicy::Threaded, typename Func>
        void execute(size_t task_count, Func&& func) const {
            APRIL_ASSERT(ptr != nullptr, "[APRIL] ExecutorRef executed before being bound to System!");
            ptr->template execute<P>(task_count, std::forward<Func>(func));
        }

    private:
        RealExecutor* ptr = nullptr;
    };
}
