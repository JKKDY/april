#pragma once

#include <chrono>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>
#include <iostream>

#include "april/monitors/monitor.hpp"

namespace april::monitor {

    struct BenchmarkResult {
        size_t steps;
        uint64_t total_updates;
        double wall_time_sec;
        double integration_time_s;
        double its_per_sec;
        double mups; // million updates per second
        double avg_step_sec;
        double median_step_sec;
        double min_step_sec;
        double max_step_sec;
        double std_dev_sec;

        void print_report() const {
            std::cout << "\n" << std::string(40, '-') << "\n";
            std::cout << " [APRIL BENCHMARK REPORT] \n";
            std::cout << std::string(40, '-') << "\n";

            std::cout << std::fixed << std::setprecision(5);
            std::cout << "  Steps processed:    " << steps << "\n";
            std::cout << "  Particles processed: " << total_updates << "\n";
            std::cout << "  Wall time (total):  " << wall_time_sec << " s\n";
            std::cout << "  Integration time:   " << integration_time_s << " s\n";
            std::cout << std::string(40, '-') << "\n";

            std::cout << std::setprecision(2);
            // Reverting to your original labels: Throughput and Performance
            std::cout << "  Throughput:         " << its_per_sec << " it/s\n";
            std::cout << "  Performance:        " << mups << " MUPS\n";
            std::cout << std::string(40, '-') << "\n";

            std::cout << std::setprecision(6);
            std::cout << "  Avg step time:      " << avg_step_sec << " s\n";
            std::cout << "  Median step time:   " << median_step_sec << " s\n";
            std::cout << "  Min step time:      " << min_step_sec << " s\n";
            std::cout << "  Max step time:      " << max_step_sec << " s\n";
            std::cout << "  Std Deviation:      " << std_dev_sec << " s\n";
            std::cout << std::string(40, '-') << "\n\n";
        }
    };

    class Benchmark : public Monitor {
    public:
       Benchmark() : Monitor(shared::Trigger::always()) {}
       explicit Benchmark(BenchmarkResult * res) : Monitor(shared::Trigger::always()), result(res) {}

       void initialize() {
          glob_start_time = std::chrono::steady_clock::now();
       }

       template<class S>
       void before_step(const core::SystemContext<S> & sys) {
          current_step_updates = sys.size();
          start_time = std::chrono::steady_clock::now();
       }

       template<class S>
       void record(const core::SystemContext<S> &) {
          end_time = std::chrono::steady_clock::now();
          const auto elapsed = std::chrono::duration<double>(end_time - start_time).count();
          timings.push_back(elapsed);
          updates += current_step_updates;
       }

       void finalize() {
          if (timings.empty()) return;

          glob_end_time = std::chrono::steady_clock::now();

          // Calculate and store the results internally
          auto res = calculate_results();
          res.print_report();

          if (result) {
             *result = res;
          }
       }

    private:
       BenchmarkResult calculate_results() {
          const double glob_total_s = std::chrono::duration<double>(glob_end_time - glob_start_time).count();
          const double total_integ_s = std::accumulate(timings.begin(), timings.end(), 0.0);
          const size_t steps = timings.size();

          std::vector<double> sorted = timings;
          std::ranges::sort(sorted);

          const double avg = total_integ_s / static_cast<double>(steps);

          double variance = 0.0;
          for (const double t : timings) variance += (t - avg) * (t - avg);

          return BenchmarkResult{
              .steps = steps,
              .total_updates = updates,
              .wall_time_sec = glob_total_s,
              .integration_time_s = total_integ_s,
              .its_per_sec = 1.0 / avg,
              .mups = (static_cast<double>(updates) / total_integ_s) / 1'000'000.0,
              .avg_step_sec = avg,
              .median_step_sec = sorted[steps / 2],
              .min_step_sec = sorted.front(),
              .max_step_sec = sorted.back(),
              .std_dev_sec = std::sqrt(variance / static_cast<double>(steps))
          };
       }

       using Clock = std::chrono::steady_clock;
       Clock::time_point glob_start_time, glob_end_time, start_time, end_time;

       std::vector<double> timings;
       uint64_t updates = 0;
       size_t current_step_updates = 0;
       BenchmarkResult * result = nullptr;
    };
}