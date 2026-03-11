#pragma once

#include <concepts>
#include <stddef.h>
namespace april::exec {

    template<typename F>
    concept IsWorkAtom = std::invocable<F, size_t>;

    template<typename T>
    concept IsExecutor = requires(T exec, size_t count) {
        // use a dummy lambda to check if the execute method exists and is templated
        { exec.execute(count, [](size_t){}) } -> std::same_as<void>;
    };

} // namespace april::exec