#pragma once

#include "april/base/backend_config.hpp"
#include "april/exec/threading/executor_concepts.hpp"
#include "april/exec/threading/executor_reference.hpp"

#if (defined(APRIL_EXECUTOR_BACKEND_OMP) +                    \
     defined(APRIL_EXECUTOR_BACKEND_NATIVE_BARRIER) +         \
     defined(APRIL_EXECUTOR_BACKEND_NATIVE_SPIN) +             \
     defined(APRIL_EXECUTOR_BACKEND_SEQUENTIAL)) > 1
    #error "[APRIL] Multiple thread executor backends selected."
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
        using DefaultThreadExecutor = SequentialExecutor;
    }

#else
    #error "[APRIL] No thread executor backend selected."
#endif

namespace april::exec {
    static_assert(IsThreadingExecutor<DefaultThreadExecutor>);
}