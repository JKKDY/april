#pragma once
#include <thread>

#ifdef __cpp_lib_hardware_interference_size
#include <new>
#endif

namespace  april::exec {
    inline const unsigned int n_threads = [] {
        const unsigned int count = std::thread::hardware_concurrency();
        return (count == 0) ? 1 : count;
    }();

#ifdef __cpp_lib_hardware_interference_size
    constexpr size_t cache_line_size = std::hardware_destructive_interference_size;
#else
    constexpr size_t cache_line_size = 64;
#endif
}