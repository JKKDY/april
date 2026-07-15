#pragma once
#include <thread>
#include <cstddef>

namespace  april::exec {
    inline const unsigned int default_thread_count = [] {
        const unsigned int count = std::thread::hardware_concurrency();
        return (count == 0) ? 1 : count;
    }();

    #ifdef __cpp_lib_hardware_interference_size
        inline constexpr size_t assumed_cache_line_size =
            std::hardware_destructive_interference_size;
    #else
        inline constexpr size_t assumed_cache_line_size = 64;
    #endif

}