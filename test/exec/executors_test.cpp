#include <gtest/gtest.h>
#include <atomic>
#include <set>
#include <mutex>

#include "april/exec/executors/context.hpp"
#include "april/exec/executors/native_barrier_executor.hpp"
#include "april/exec/executors/native_spin_executor.hpp"
#include "april/exec/executors/omp_executor.hpp"

using namespace april;
using namespace april::exec;

template <typename Executor>
class ExecutorTest : public testing::Test {
protected:
    Executor::Config config;
    // We'll use 4 threads for testing to ensure we actually trigger parallelism
    ExecutorTest() { config.n_threads = 4; }
};

// Define the list of executors to test
#ifdef _OPENMP
using ExecutorTypes = testing::Types<NativeBarrierExecutor, NativeSpinExecutor, OmpExecutor>;
#else
using ExecutorTypes = testing::Types<NativeBarrierExecutor, NativeSpinExecutor>;
#endif
TYPED_TEST_SUITE(ExecutorTest, ExecutorTypes);

// 1. Test Work Completion & Index Range
TYPED_TEST(ExecutorTest, ExecutesAllTasksWithCorrectIndices) {
    TypeParam executor(this->config);

    const size_t num_threads = executor.num_threads();
    const size_t num_tasks = num_threads * 4;

    std::atomic<size_t> counter{0};
    std::mutex set_mutex;
    std::set<int> unique_thread_ids;

    executor.template execute<april::ParallelPolicy::Threaded>(num_tasks, [&](size_t ) {
        int tid = thread_index();

        std::this_thread::sleep_for(std::chrono::milliseconds(1));
        
        // Check ID range
        EXPECT_GE(tid, 0);
        EXPECT_LT(tid, static_cast<int>(executor.num_threads()));
        
        {
            std::lock_guard lock(set_mutex);
            unique_thread_ids.insert(tid);
        }
        counter.fetch_add(1);
    });

    EXPECT_EQ(counter.load(), num_tasks);
    // Ensure we actually used more than 1 thread
    EXPECT_GT(unique_thread_ids.size(), 1u);
}

// 2. Test Serial Policy (thread_index should be -1 and is_parallel should be false)
TYPED_TEST(ExecutorTest, SerialPolicyDoesNotSetContext) {
    TypeParam executor(this->config);
    
    // Outside execute
    EXPECT_FALSE(is_parallel());
    EXPECT_EQ(thread_index(), -1);

    executor.template execute<ParallelPolicy::Serial>(10, [&](size_t) {
        // Inside Serial execute (User's specific request: no context in serial)
        EXPECT_FALSE(is_parallel());
        EXPECT_EQ(thread_index(), -1);
    });

    // Back outside
    EXPECT_FALSE(is_parallel());
}

// 3. Test RAII Nesting Logic
TYPED_TEST(ExecutorTest, HandlesNestedContextsCorrectly) {
    TypeParam executor(this->config);

    // Initial state
    EXPECT_EQ(thread_index(), -1);

    executor.execute(1, [&](size_t ) {
        int outer_tid = thread_index();
        EXPECT_GE(outer_tid, 0);

        {
            // Simulate a nested region or manual guard
            internal::ScopedThreadContext nested_guard(99);
            EXPECT_EQ(thread_index(), 99);
        }

        // RAII must restore the ID assigned by the executor
        EXPECT_EQ(thread_index(), outer_tid);
    });

    EXPECT_EQ(thread_index(), -1);
}

// 4. Test Zero Batch Count (Should return immediately)
TYPED_TEST(ExecutorTest, HandlesZeroBatchCount) {
    TypeParam executor(this->config);
    bool called = false;
    executor.execute(0, [&](size_t) { called = true; });
    EXPECT_FALSE(called);
}