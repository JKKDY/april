#pragma once


#define AP_EXECUTOR_BACKEND_OMP



#if defined(AP_EXECUTOR_BACKEND_OMP)
#include "april/exec/executors/omp_executor.hpp"

namespace april::exec {
    using Executor = OmpExecutor;
}

#elif defined(AP_EXECUTOR_BACKEND_NATIVE)
namespace april::exec {
    using Executor = NativeExecutor;
}
#elif defined(AP_EXECUTOR_BACKEND_SEQUENTIAL)
#include "executors/sequential_executor.hpp"
namespace april::exec {
    using Executor = ::april::exec::SequentialExecutor;
}
#endif


static_assert(april::exec::IsExecutor<april::exec::Executor>);


