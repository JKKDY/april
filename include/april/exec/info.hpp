#pragma once
#include <thread>

namespace  april::exec {
    inline const unsigned int N_CPU_THREADS = [] {
        const unsigned int count = std::thread::hardware_concurrency();
        return (count == 0) ? 1 : count;
    }();

    constexpr size_t CACHE_LINE_SIZE = 64;
}