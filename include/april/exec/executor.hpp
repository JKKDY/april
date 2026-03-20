#pragma once


#define AP_EXECUTOR_USE_NATIVE_SPIN

#if defined(AP_EXECUTOR_USE_OMP)
#include "april/exec/executors/omp_executor.hpp"

namespace april::exec {
    using Executor = OmpExecutor;
}

#elif defined(AP_EXECUTOR_USE_NATIVE_BARRIER)

#include "april/exec/executors/native_executor.hpp"

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


