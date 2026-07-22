#pragma once
#include <thread>
#include <cstddef>

namespace  april::exec {
    inline const unsigned int default_thread_count = [] {
        const unsigned int count = std::thread::hardware_concurrency();
        return (count == 0) ? 1 : count;
    }();

    inline constexpr size_t assumed_cache_line_size = 64;

}