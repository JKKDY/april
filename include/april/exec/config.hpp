#pragma once

#include "april/base/config.hpp"
#include "april/exec/policy.hpp"


#if defined(AP_EXECUTOR_USE_OMP)
#include "april/exec/executors/omp_executor.hpp"

namespace april::exec {
    using Executor = OmpExecutor;
}

#elif defined(AP_EXECUTOR_USE_NATIVE_BARRIER)

#include "april/exec/executors/native_barrier_executor.hpp"

namespace april::exec {
    using Executor = NativeBarrierExecutor;
}
#elif defined(AP_EXECUTOR_USE_NATIVE_SPIN)

#include "april/exec/executors/native_spin_executor.hpp"

namespace april::exec {
    using Executor = NativeSpinExecutor;
}

#elif defined(AP_EXECUTOR_BACKEND_SEQUENTIAL)

#include "executors/sequential_executor.hpp"
namespace april::exec {
    using Executor = ::april::exec::SequentialExecutor;
}
#endif


static_assert(april::exec::IsExecutor<april::exec::Executor>);

namespace april {
    template <ParallelPolicy P = ParallelPolicy::Threaded, VectorPolicy V =VectorPolicy::Auto>
    struct ExecutionConfig {
        static constexpr VectorPolicy vector_policy = V;
        static constexpr ParallelPolicy parallel_policy = P;
        exec::Executor thread_executor = exec::Executor(exec::N_CPU_THREADS);
    };

    namespace exec {
        template <typename C>
        concept IsExecutionConfig =
         requires (C c) {
            { C::vector_policy } -> std::convertible_to<VectorPolicy>;
            { C::parallel_policy } -> std::convertible_to<ParallelPolicy>;
            { c.thread_executor } -> IsExecutor;
         };
    }
}


