#include <gtest/gtest.h>
#include <set>
#include <vector>
#include <unordered_set>

#include "april/containers/batching/topology_batch.hpp"
#include "april/exec/parallel_utils.hpp"
#include "april/math/range.hpp"

// Dummy config to force deterministic testing
struct TestConfig {
    size_t fixed_blocks;
    size_t alignment;

    size_t calculate_num_blocks(size_t /* elems */) const {
        return fixed_blocks;
    }
};

// Helper for Google Test to print ranges nicely if a test fails
namespace april::math {
    inline bool operator==(const Range& a, const Range& b) {
        return a.start == b.start && a.stop == b.stop;
    }
}


TEST(LinearSchedule, PerfectDivision) {
    const auto blocks = april::exec::make_linear_schedule(april::math::Range{0, 100}, TestConfig{4, 1});

    ASSERT_EQ(blocks.size(), 4);
    EXPECT_EQ(blocks[0], (april::math::Range{0, 25}));
    EXPECT_EQ(blocks[1], (april::math::Range{25, 50}));
    EXPECT_EQ(blocks[2], (april::math::Range{50, 75}));
    EXPECT_EQ(blocks[3], (april::math::Range{75, 100}));
}

TEST(LinearSchedule, RemainderHandling) {
    // 10 elements into 3 blocks. Should distribute: 4, 3, 3
    const auto blocks = april::exec::make_linear_schedule(april::math::Range{0, 10}, TestConfig{3, 1});

    ASSERT_EQ(blocks.size(), 3);
    EXPECT_EQ(blocks[0], (april::math::Range{0, 4}));
    EXPECT_EQ(blocks[1], (april::math::Range{4, 7}));
    EXPECT_EQ(blocks[2], (april::math::Range{7, 10}));
}

TEST(LinearSchedule, AlignmentPacking) {
    // 10 elements, 2 blocks, but aligned to 4 (e.g., AVX width).
    // Total packed = 2 (8 elements). Vectors per block = 1.
    // Block 0 gets 1 vector (4 elements).
    // Block 1 gets 1 vector (4 elements) + scalar tail (2 elements).
    const auto blocks = april::exec::make_linear_schedule(april::math::Range{0, 10}, TestConfig{2, 4});

    ASSERT_EQ(blocks.size(), 2);
    EXPECT_EQ(blocks[0], (april::math::Range{0, 4}));
    EXPECT_EQ(blocks[1], (april::math::Range{4, 10}));
}

TEST(SymmetricSchedule, CoverageAndIndependence) {
    // 4 blocks = 4 diagonals (Phase 0), and 3 off-diagonal phases.
    auto schedule = april::exec::make_symmetric_schedule(april::math::Range{0, 100}, TestConfig{4, 1});

    // 1. Check Phase counts
    EXPECT_EQ(schedule.diagonals.size(), 4);
    EXPECT_EQ(schedule.off_diagonals.size(), 3);

    std::set<std::pair<size_t, size_t>> seen_pairs;

    // Track diagonals
    for (const auto& diag : schedule.diagonals) {
        seen_pairs.insert({diag.start, diag.start});
    }

    // 2. Check Independence inside each off-diagonal phase
    for (const auto& phase : schedule.off_diagonals) {
        std::set<size_t> active_blocks_in_phase;
        for (const auto& [b1, b2] : phase) {
            // Ensure no block is touched twice in the same phase (Race Condition Check!)
            EXPECT_TRUE(active_blocks_in_phase.insert(b1.start).second);
            EXPECT_TRUE(active_blocks_in_phase.insert(b2.start).second);

            seen_pairs.insert({b1.start, b2.start});
        }
    }

    // 3. Check Total Coverage (Gauss Sum)
    // 4 blocks -> 4 + 3 + 2 + 1 = 10 unique pairs
    EXPECT_EQ(seen_pairs.size(), 10);
}

TEST(SymmetricSchedule, OddBlockCountCorrection) {
    // 3 blocks requested -> should internally pad to 4 (even requirement)
    auto schedule = april::exec::make_symmetric_schedule(april::math::Range{0, 30}, TestConfig{3, 1});

    // It should have 3 valid diagonals (the 4th is empty and filtered out)
    EXPECT_EQ(schedule.diagonals.size(), 3);
    // It should still have 3 phases (padded to 4 blocks: 4-1 = 3)
    EXPECT_EQ(schedule.off_diagonals.size(), 3);
}

TEST(BipartiteSchedule, CoverageAndIndependence) {
    const auto schedule = april::exec::make_bipartite_schedule(
        april::math::Range{0, 100},  // Range 1
        april::math::Range{100, 200}, // Range 2
        TestConfig{4, 1}
    );

    // B blocks -> exactly B phases
    ASSERT_EQ(schedule.phases.size(), 4);

    std::set<std::pair<size_t, size_t>> seen_pairs;

    for (const auto& phase : schedule.phases) {
        std::set<size_t> active_b1, active_b2;

        for (const auto& [b1, b2] : phase) {
            // No block from Range 1 or 2 should be used more than once per phase
            EXPECT_TRUE(active_b1.insert(b1.start).second);
            EXPECT_TRUE(active_b2.insert(b2.start).second);

            seen_pairs.insert({b1.start, b2.start});
        }
    }

    // Total coverage: 4 blocks in R1 * 4 blocks in R2 = 16 unique interaction pairs
    EXPECT_EQ(seen_pairs.size(), 16);
}


TEST(ScheduleEdgeCases, EmptyRanges) {
    TestConfig config{4, 1};

    auto linear = april::exec::make_linear_schedule(april::math::Range{0, 0}, config);
    EXPECT_EQ(linear.size(), 1); // Fallback safety
    EXPECT_EQ(linear[0], (april::math::Range{0, 0}));

    auto symmetric = april::exec::make_symmetric_schedule(april::math::Range{0, 0}, config);
    EXPECT_TRUE(symmetric.diagonals.empty());
    EXPECT_TRUE(symmetric.off_diagonals.empty());

    auto bipartite = april::exec::make_bipartite_schedule(
        april::math::Range{0, 0},
        april::math::Range{100, 100},
        config
    );
    EXPECT_TRUE(bipartite.phases.empty());
}

TEST(ScheduleEdgeCases, TinyRange) {
    // 2 elements, but forcing 4 blocks.
    // The algorithm front-loads remainders, so blocks 0 and 1 will get the work,
    // and blocks 2 and 3 will be empty.
    const auto blocks = april::exec::make_linear_schedule(april::math::Range{0, 2}, TestConfig{4, 1});

    ASSERT_EQ(blocks.size(), 4);
    EXPECT_EQ(blocks[0], (april::math::Range{0, 1}));
    EXPECT_EQ(blocks[1], (april::math::Range{1, 2}));
    EXPECT_EQ(blocks[2], (april::math::Range{2, 2}));
    EXPECT_EQ(blocks[3], (april::math::Range{2, 2}));
}

TEST(ScheduleEdgeCases, ImbalancedBipartite) {
    // Range 1 is massive (100 elements), Range 2 is tiny (2 elements)
    auto schedule = april::exec::make_bipartite_schedule(
        april::math::Range{0, 100},
        april::math::Range{100, 102},
        TestConfig{4, 1}
    );

    // Because Range 2 only has 2 elements but is split into 4 blocks,
    // 2 of those blocks will be empty [start, start].
    // Your bipartite scheduler has a filter: `if (blocks[i].start < blocks[i].stop)`
    // We must verify that no empty blocks made it into the final phases.

    for (const auto& phase : schedule.phases) {
        for (const auto& [b1, b2] : phase) {
            EXPECT_LT(b1.start, b1.stop);
            EXPECT_LT(b2.start, b2.stop);
        }
    }
}



using namespace april::container::batching;
using Node = april::ParticleID;
using Pairs = april::utility::graph::EdgeList<Node>;

// Helper: Count total pairs in the output schedule
size_t count_total_pairs(const std::vector<std::vector<Pairs>>& phases) {
    size_t count = 0;
    for (const auto& phase : phases) {
        for (const auto& batch : phase) {
            count += batch.size();
        }
    }
    return count;
}

// Helper: Assert that a phase has absolutely zero overlapping nodes between its batches
void assert_phase_independence(const std::vector<Pairs>& phase) {
    std::unordered_set<Node> seen_nodes_in_phase;

    for (const auto& batch : phase) {
        std::unordered_set<Node> nodes_in_batch;
        // Collect all unique nodes used by this specific batch
        for (const auto& [u, v] : batch) {
            nodes_in_batch.insert(u);
            nodes_in_batch.insert(v);
        }

        // Ensure none of these nodes are being touched by another batch in this phase
        for (const Node n : nodes_in_batch) {
            ASSERT_TRUE(seen_nodes_in_phase.insert(n).second)
                << "Cross-batch race condition detected on Node " << n;
        }
    }
}


TEST(TopologyScheduling, EmptyInput) {
    std::vector<Pairs> input = {};
    auto phases = build_concurrent_phases(input, 100, 1);
    EXPECT_TRUE(phases.empty());

    // Vector with empty edge lists
    std::vector<Pairs> empty_lists = {{}, {}};
    phases = build_concurrent_phases(empty_lists, 100, 1);
    EXPECT_TRUE(phases.empty());
}

TEST(TopologyScheduling, PerfectParallelism) {
    // 3 completely independent interactions
    const std::vector<Pairs> input = {
        {{0, 1}},
        {{2, 3}},
        {{4, 5}}
    };

    // max_partition_size = 10, min_threshold = 1
    const auto phases = build_concurrent_phases(input, 10, 1);

    // Because they are independent, they should all be packed into a single Phase
    ASSERT_EQ(phases.size(), 1);
    // There should be 3 distinct atomic batches in this phase
    EXPECT_EQ(phases[0].size(), 3);
    EXPECT_EQ(count_total_pairs(phases), 3);
    assert_phase_independence(phases[0]);
}

TEST(TopologyScheduling, ConnectedComponentGrouping) {
    // A "Hub and Spoke" molecule inside a single topology group
    // Node 0 is connected to 1, 2, and 3.
    const std::vector<Pairs> input = {
        {{0, 1}, {0, 2}, {0, 3}, {4, 5}}
    };

    const auto phases = build_concurrent_phases(input, 10, 1);

    ASSERT_EQ(phases.size(), 1);
    // The disjoint-set-union should group {0,1,2,3} into ONE atomic batch,
    // and {4,5} into a SECOND atomic batch.
    ASSERT_EQ(phases[0].size(), 2);
    EXPECT_EQ(count_total_pairs(phases), 4);
    assert_phase_independence(phases[0]);
}

TEST(TopologyScheduling, ConflictResolution) {
    // Two separate interaction types (e.g., Bonds and Angles)
    // They share node '2', creating a cross-topology conflict.
    const std::vector<Pairs> input = {
        {{0, 1}, {1, 2}}, // Topology A
        {{2, 3}, {4, 5}}  // Topology B
    };

    const auto phases = build_concurrent_phases(input, 10, 1);

    // Because Node 2 is contested between Top A and Top B, they CANNOT be in the same phase.
    // However, the isolated pair {4, 5} from Top B should pack cleanly with Top A's batch.
    ASSERT_EQ(phases.size(), 2);
    EXPECT_EQ(count_total_pairs(phases), 4);

    for (const auto& phase : phases) {
        assert_phase_independence(phase);
    }
}


TEST(TopologyScheduling, SequentialFallbackFlattening) {
    // 5 conflicting pairs -> The graph will be forced to create 5 separate phases
    // because they all share node '0'.
    const std::vector<Pairs> input = {
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{0, 4}}, {{0, 5}}
    };

    // min_threshold = 3.
    // The greedy colorer will make 5 phases, each with exactly 1 batch.
    // ALL of them are below the threshold of 3.
    // Therefore, ALL of them should be dumped into the sequential fallback.
    const auto phases = build_concurrent_phases(input, 10, 3);

    // There should be exactly ONE phase outputted (the sequential fallback)
    ASSERT_EQ(phases.size(), 1);

    // Inside that phase, there must be exactly ONE flattened batch
    ASSERT_EQ(phases[0].size(), 1);

    // That single batch must contain all 5 pairs
    EXPECT_EQ(phases[0][0].size(), 5);
    EXPECT_EQ(count_total_pairs(phases), 5);

    // NOTE: We DO NOT call assert_phase_independence(phases[0]) here!
    // The sequential fallback is inherently conflicted, which is mathematically
    // safe because it is dispatched to a single thread.
}

TEST(TopologyScheduling, MixedParallelAndSequential) {
    // 3 independent pairs.
    // 3 conflicting pairs (all share Node 0).
    const std::vector<Pairs> input = {
        {{10, 11}}, {{12, 13}}, {{14, 15}},
        {{0, 1}}, {{0, 2}}, {{0, 3}}
    };

    // min_threshold = 2
    const auto phases = build_concurrent_phases(input, 10, 2);

    // HOW THE ENGINE ACTUALLY PACKS THIS:
    // Phase 0: Takes {0,1} AND all 3 independent batches. Total = 4 batches. (Parallel)
    // Phase 1: Takes {0,2}. Total = 1 batch. (Fails threshold)
    // Phase 2: Takes {0,3}. Total = 1 batch. (Fails threshold)
    // Fallback: Merges Phase 1 and Phase 2 into a single sequential batch.

    ASSERT_EQ(phases.size(), 2);

    // Verify Phase 0 (Parallel - 4 distinct batches)
    EXPECT_EQ(phases[0].size(), 4);
    assert_phase_independence(phases[0]);

    // Verify Phase 1 (Sequential Fallback)
    EXPECT_EQ(phases[1].size(), 1); // Flattened to exactly 1 batch
    EXPECT_EQ(phases[1][0].size(), 2); // Containing the 2 leftover conflicting pairs ({0,2} and {0,3})
    EXPECT_EQ(count_total_pairs(phases), 6);
}