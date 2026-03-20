#include "april/exec/executors/native_executor.hpp"
#include "april/exec/executors/omp_executor.hpp"
#include "../include/april/exec/executors/performance_executor.hpp"

#include <iostream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <string>
#include <iomanip>

using namespace april::exec;

template <typename Executor, typename Task>
void run_benchmark(Executor& exec, Task&& task, const std::string& label, int N = 10000) {
    std::vector<double> timings;
    timings.reserve(N);

    // Warm-up: discard the first few runs to prime caches/JIT
    for (int i = 0; i < 10; ++i) {
        exec.execute(1, task);
    }

    for (int i = 0; i < N; i++) {
        const auto start = std::chrono::steady_clock::now();

        exec.execute(1, task);

        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration<double> elapsed = end - start;
        timings.push_back(elapsed.count());
    }

    // --- Statistics ---
    std::ranges::sort(timings);
    const double sum = std::accumulate(timings.begin(), timings.end(), 0.0);
    const double avg = sum / N;
    const double median = timings[N / 2];
    const double p95 = timings[static_cast<size_t>(N * 0.95)];

    const double sq_sum = std::inner_product(timings.begin(), timings.end(), timings.begin(), 0.0);
    const double stdev = std::sqrt(std::abs(sq_sum / N - avg * avg));

    // --- Output Table Format ---
    std::cout << "\n=== " << label << " (" << N << " iterations) ===\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Average: " << avg << "s\n";
    std::cout << "Median:  " << median << "s\n";
    std::cout << "StdDev:  " << stdev << "s\n";
    std::cout << "P95:     " << p95 << "s\n";
    std::cout << "Min/Max: " << timings.front() << "s / " << timings.back() << "s\n";
}

int main() {
    auto task = [](size_t i) {
        double x = static_cast<double>(i);
        x = x * 1.0000001 + 0.0000001;
        asm volatile("" : "+x"(x));
    };

    constexpr int iterations = 10000;

    NativeExecutor n;
    run_benchmark(n, task, "Native Executor", iterations);

    OmpExecutor o;
    run_benchmark(o, task, "OpenMP Executor", iterations);

    return 0;
}