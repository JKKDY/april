#pragma once

#include <concepts>

#include "april/exec/policy.hpp"
#include "april/exec/threading/executor_config.hpp"
#include "april/exec/threading/executor_concepts.hpp"




// public facing configuration helpers
namespace april {

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