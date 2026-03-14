#pragma once

#include <vector>
#include <algorithm>
#include <ranges>
#include <numeric>

#include "april/base/types.hpp"
#include "april/math/range.hpp"
#include "april/exec/info.hpp"

namespace april::exec {

    template<typename C>
    concept IsBlockConfig = requires (const C c, size_t elems){
        {c.calculate_blocks(elems)} -> std::convertible_to<size_t>;
        {c.alignment} -> std::convertible_to<size_t>;
    };

    struct BlockConfig {
        size_t min_tasks;
        size_t target_elements_per_task;
        size_t alignment;

        BlockConfig() : min_tasks(1), target_elements_per_task(256), alignment(1) {}

        explicit BlockConfig(const size_t n_threads,  const size_t oversubscription = 4, const size_t target_size = 256)
            : min_tasks(n_threads * oversubscription),
              target_elements_per_task(target_size),
              alignment(CACHE_LINE_SIZE / sizeof(vec3::type)) {}  // aligned to cache line width

        // heuristic for determining optimal partition count
        [[nodiscard]] size_t calculate_blocks(const size_t total_elements) const {
            if (total_elements == 0) return 0;

            const size_t ideal_B = total_elements / std::max<size_t>(1, target_elements_per_task);

            size_t B = std::max(ideal_B, min_tasks);

            const size_t total_packed = total_elements / std::max<size_t>(1, alignment);
            if (total_packed < B) {
                B = std::max<size_t>(1, total_packed);
            }

            return B;
        }
    };

    struct SymmetricPhaseSchedule {
        std::vector<math::Range> diagonals; // Phase 0
        std::vector<std::vector<std::pair<math::Range, math::Range>>> off_diagonals; // Phases 1 to B-1
    };

    struct BipartitePhaseSchedule {
        std::vector<std::vector<std::pair<math::Range, math::Range>>> phases;
    };


    // partition a linear index range into an exact number of blocks
    [[nodiscard]] inline std::vector<math::Range> make_linear_schedule(
        const math::Range& range,
        const size_t B,
        const size_t alignment)
    {
        if (range.empty() || B == 0) {
            return std::vector<math::Range>(std::max<size_t>(1, B), {range.start, range.start});
        }

        const size_t total_elements = range.size();
        const size_t total_packed = total_elements / alignment;

        const size_t vectors_per_block = total_packed / B;
        const size_t remainder_vectors = total_packed % B;

        std::vector<math::Range> blocks(B);
        size_t current = range.start;

        for (size_t i = 0; i < B; ++i) {
            const size_t v_count = vectors_per_block + (i < remainder_vectors ? 1 : 0);
            size_t size = v_count * alignment;

            // Attach scalar tail to the very last block
            if (i == B - 1) {
                size += total_elements % alignment;
            }

            blocks[i] = {current, current + size};
            current += size;
        }

        return blocks;
    }

    // partition a linear index range using block configuration heuristics
    [[nodiscard]] std::vector<math::Range> make_linear_schedule(
        const math::Range& range,
        IsBlockConfig auto const& config)
    {
        const size_t total_elements = range.size();
        if (total_elements == 0) {
            return {{range.start, range.start}};
        }

        const size_t B = config.calculate_blocks(total_elements);

        if (B <= 1) {
            return {range};
        }

        return make_linear_schedule(range, B, config.alignment);
    }




    // Schedule pairwise self-interactions into independent execution phases
    [[nodiscard]] SymmetricPhaseSchedule make_symmetric_schedule(
        const math::Range& range,
        IsBlockConfig auto const& config)
    {
        SymmetricPhaseSchedule schedule;
        if (range.empty()) return schedule;

        // partition range
        std::vector<math::Range> blocks = make_linear_schedule(range, config);
        size_t B = blocks.size();

        // make B even (requirement for round-robin tournament)
        if (B % 2 != 0) {
            blocks.push_back({range.stop, range.stop});
            ++B;
        }

        // Phase 0: Diagonals
        schedule.diagonals.reserve(B);
        for (size_t i = 0; i < B; ++i) {
            if (blocks[i].start < blocks[i].stop) {
                schedule.diagonals.push_back(blocks[i]);
            }
        }

        // Phases 1 to B-1: Round-Robin Tournament
        schedule.off_diagonals.resize(B - 1);
        std::vector<size_t> circle(B);
        std::iota(circle.begin(), circle.end(), 0);

        for (size_t phase = 0; phase < B - 1; ++phase) {
            schedule.off_diagonals[phase].reserve(B / 2);
            for (size_t i = 0; i < B / 2; ++i) {
                size_t b1 = circle[i];
                size_t b2 = circle[B - 1 - i];

                const size_t idx1 = std::min(b1, b2);
                const size_t idx2 = std::max(b1, b2);

                if (blocks[idx1].start < blocks[idx1].stop &&
                    blocks[idx2].start < blocks[idx2].stop)
                {
                    schedule.off_diagonals[phase].emplace_back(blocks[idx1], blocks[idx2]);
                }
            }

            // Rotate circle array: keep index 0 fixed, rotate others right
            const size_t last = circle.back();
            for (size_t i = B - 1; i > 1; --i) circle[i] = circle[i - 1];
            circle[1] = last;
        }

        return schedule;
    }

    // schedule pair interactions between two bipartite ranges into independent execution phases
    [[nodiscard]] BipartitePhaseSchedule make_bipartite_schedule(
        const math::Range& range1,
        const math::Range& range2,
        IsBlockConfig auto const& config)
    {
        BipartitePhaseSchedule schedule;
        if (range1.empty() || range2.empty()) return schedule;

        // determine common block count B based on the larger range
        const size_t max_elements = std::max(range1.size(), range2.size());
        size_t B = config.calculate_blocks(max_elements);

        if (B == 0) B = 1; // Fallback safety

        // partition both ranges into exactly B blocks
        std::vector<math::Range> blocks1 = make_linear_schedule(range1, B, config.alignment);
        std::vector<math::Range> blocks2 = make_linear_schedule(range2, B, config.alignment);

        AP_ASSERT(blocks1.size() == blocks2.size(), "Bipartite blocks mismatch");

        // cyclic Diagonal Scheduling (Phase p: Pair block i with block (i+p)%B)
        schedule.phases.resize(B);
        for (size_t p = 0; p < B; ++p) {
            schedule.phases[p].reserve(B);
            for (size_t i = 0; i < B; ++i) {
                const size_t j = (i + p) % B;

                if (blocks1[i].start < blocks1[i].stop &&
                    blocks2[j].start < blocks2[j].stop)
                {
                    schedule.phases[p].emplace_back(blocks1[i], blocks2[j]);
                }
            }
        }

        return schedule;
    }

} // namespace april::exec