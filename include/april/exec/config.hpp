#pragma once

#include <concepts>

#include "april/base/backend_config.hpp"
#include "april/exec/policy.hpp"
#include "april/exec/thread_executor.hpp"
#include "april/exec/threading/executor_concepts.hpp"


// threading executor backend selection
#if (defined(APRIL_EXECUTOR_BACKEND_OMP) + \
defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER) + \
defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN) + \
defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)) > 1
#error "[APRIL] Multiple executor backends defined. Select exactly one."
#endif


#if defined(APRIL_EXECUTOR_BACKEND_OMP)
    #include "april/exec/threading/backends/omp_executor.hpp"
    namespace april::exec {
        using DefaultThreadExecutor = OmpExecutor;
    }
#elif defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER)
#include "april/exec/threading/backends/native_barrier_executor.hpp"
namespace april::exec {
        using DefaultThreadExecutor = NativeBarrierExecutor;
    }
#elif defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN)
#include "april/exec/threading/backends/native_spin_executor.hpp"
namespace april::exec {
        using DefaultThreadExecutor = NativeSpinExecutor;
    }
#elif defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)
#include "april/exec/threading/backends/sequential_executor.hpp"
namespace april::exec {
        using DefaultThreadExecutor = ::april::exec::SequentialExecutor;
    }
#else
#error "[APRIL] No executor backend selected"
#endif




// public facing configuration helpers
namespace april {
    static_assert(april::exec::IsThreadingExecutor<exec::DefaultThreadExecutor>);


    /**
    * @brief Rebindable, non-owning reference to a thread executor.
    */
    template<exec::IsThreadingExecutor E = exec::DefaultThreadExecutor>
    struct RuntimeConfig {
        using ThreadExecutor = E;
        using ExecutorConfig = E::Config;

        ExecutorConfig executor_config{};
    };

    template<
        ParallelPolicy Parallel = ParallelPolicy::Threaded,
        VectorPolicy Vector = VectorPolicy::Auto
    >
    struct CompileTimeConfig {
        static constexpr VectorPolicy vector_policy = Vector;
        static constexpr ParallelPolicy parallel_policy = Parallel;
    };

    struct ExecutionConfig : RuntimeConfig<>, CompileTimeConfig<> {};

} // namespace april


namespace april::exec {

    template<typename C>
    concept IsExecutionConfig =
        requires(C& config) {
        typename C::ThreadExecutor;
        typename C::ExecutorConfig;

        requires IsThreadingExecutor<typename C::ThreadExecutor>;

        { C::vector_policy } -> std::convertible_to<VectorPolicy>;
        { C::parallel_policy } -> std::convertible_to<ParallelPolicy>;
        { config.executor_config };
    };
} // namespace april::exec