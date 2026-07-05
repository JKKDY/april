#pragma once

#include <concepts>
#include <stddef.h>

namespace april::exec {
    template<typename F>
    concept IsWorkAtom = std::invocable<F, size_t>;

    template<typename T>
    concept IsExecutor = requires(const T exec, size_t count) {

        // execute takes in the number of batches and a callable with the current index of the batch being processed
        { exec.execute(count, [](size_t){}) } -> std::same_as<void>;

        // return number of threads of the executor
        { exec.num_threads() } -> std::convertible_to<size_t>;
    };
} // namespace april::exec