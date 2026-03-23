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

    template<exec::IsExecutor E = exec::Executor>
    struct RunTimeConfig {
        using ThreadExecutor = E;
        using ExecutorConfig = E::Config;
        ExecutorConfig executer_config;
        // size_t task_chunk_size = 256; // future stub
    };

    enum class Policy {
        Enabled,
        Disabled,
        Auto
    };

    template<
        ParallelPolicy P = ParallelPolicy::Threaded,
        VectorPolicy V = VectorPolicy::Auto>
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
    }
}


