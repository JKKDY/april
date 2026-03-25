#pragma once

#include "april/base/config.hpp"
#include "april/exec/policy.hpp"
#include "executors/executor_traits.hpp"



#if (defined(APRIL_EXECUTOR_BACKEND_OMP) + \
defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER) + \
defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN) + \
defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)) > 1

#error "[APRIL] Multiple executor backends defined. Select exactly one."

#endif


#if defined(APRIL_EXECUTOR_BACKEND_OMP)

    #include "april/exec/executors/omp_executor.hpp"
    namespace april::exec {
        using Executor = OmpExecutor;
    }

#elif defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER)

#include "april/exec/executors/native_barrier_executor.hpp"
namespace april::exec {
        using Executor = NativeBarrierExecutor;
    }

#elif defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN)

#include "april/exec/executors/native_spin_executor.hpp"
namespace april::exec {
        using Executor = NativeSpinExecutor;
    }

#elif defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)

#include "executors/sequential_executor.hpp"
namespace april::exec {
        using Executor = ::april::exec::SequentialExecutor;
    }

#else
#error "[APRIL] No executor backend selected"
#endif


static_assert(april::exec::IsExecutor<april::exec::Executor>);

namespace april {

    template<exec::IsExecutor E = exec::Executor>
    struct RunTimeConfig {
        using ThreadExecutor = E;
        using ExecutorConfig = E::Config;
        ExecutorConfig executer_config;
    };

    template<ParallelPolicy P = ParallelPolicy::Threaded, VectorPolicy V = VectorPolicy::Auto>
    struct CompileTimeConfig {
        static constexpr VectorPolicy vector_policy = V;
        static constexpr ParallelPolicy parallel_policy = P;
    };

    struct ExecutionConfig : RunTimeConfig<>, CompileTimeConfig<> {};

    namespace exec {
        template <typename C>
        concept IsExecutionConfig = requires (C c) {
            { C::vector_policy } -> std::convertible_to<VectorPolicy>;
            { C::parallel_policy } -> std::convertible_to<ParallelPolicy>;
            { c.executer_config };
        } && IsExecutor<typename C::ThreadExecutor>;


        template <typename RealExecutor>
        class ExecutorRef {
        public:
            using Config =  RealExecutor::Config;

            ExecutorRef() = default;
            explicit ExecutorRef(RealExecutor* p) : ptr(p) {}

            void bind(RealExecutor* p) {
                ptr = p;
            }

            [[nodiscard]] size_t num_threads() const noexcept {
                AP_ASSERT(ptr != nullptr, "[APRIL] ExecutorRef executed before being bound to System!");
                return ptr->num_threads();
            }

            template <ParallelPolicy P = ParallelPolicy::Threaded, typename Func>
            void execute(size_t task_count, Func&& func) const {
                AP_ASSERT(ptr != nullptr, "[APRIL] ExecutorRef executed before being bound to System!");
                ptr->template execute<P>(task_count, std::forward<Func>(func));
            }

        private:
            RealExecutor* ptr = nullptr;
        };
    }
}


