#pragma once
#include <thread>

namespace  april::exec {
    inline const unsigned int CPU_THREADS = [] {
        const unsigned int count = std::thread::hardware_concurrency();
        return (count == 0) ? 1 : count;
    }();
}