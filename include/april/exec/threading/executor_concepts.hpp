#pragma once

#include <concepts>

#include "april/exec/policy.hpp"

namespace april::exec {
    template<typename F>
    concept IsIndexedWork = std::invocable<F&, std::size_t>;

    template<typename T>
    concept IsThreadingExecutor =
    requires {
        typename T::Config;
    } &&
    requires(
        const T& executor,
        const typename T::Config& config,
        std::size_t count
    ) {
        T{config};

        // execute takes in the number of batches and a callable with the current index of the batch being processed
        {
            executor.template execute<ParallelPolicy::Serial>(
                count, [](std::size_t) {}
            )
        } -> std::same_as<void>;

        {
            executor.template execute<ParallelPolicy::Threaded>(
                count, [](std::size_t) {}
            )
        } -> std::same_as<void>;


        // return number of threads of the executor
        { executor.num_threads() } -> std::convertible_to<std::size_t>;
    };

} // namespace april::exec