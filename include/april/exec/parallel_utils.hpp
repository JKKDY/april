#pragma once

#include <vector>
#include <algorithm>
#include <span>
#include <ranges>
#include <numeric>

#include "april/base/types.hpp"
#include "april/math/range.hpp"


namespace april::exec {

    struct BlockConfig {
        BlockConfig() = default;
        explicit BlockConfig(const unsigned n_threads, const unsigned oversubscription = 4):
        target_chunk_size(n_threads * oversubscription) {}

        size_t target_chunk_size = 256;
        size_t alignment = packed::size();


    };

    struct SymmetricPhaseSchedule {
        std::vector<math::Range> diagonals; // Phase 0
        std::vector<std::vector<std::pair<math::Range, math::Range>>> off_diagonals; // Phases 1 to B-1
    };

    struct BipartitePhaseSchedule {
        std::vector<std::vector<std::pair<math::Range, math::Range>>> phases;
    };



    // partition a linear index range into equally sized work blocks
    [[nodiscard]] inline std::vector<math::Range> make_linear_schedule (
        const math::Range& range,
        const BlockConfig& config = {},
        std::optional<size_t> force_B = std::nullopt)
    {
        const size_t total_elements = range.size();
        const size_t v_size = config.alignment;
        const size_t total_packed = total_elements / v_size;

        // 1. Determine target block count
        size_t B = force_B.value_or(std::max<size_t>(1, config.target_chunk_size));

        // Clamp chunks only if we aren't forced (prevents out-of-bounds in bipartite math)
        if (!force_B.has_value()) {
            if (total_packed < B) {
                B = std::max<size_t>(1, total_packed);
            }
        }

        // 2. Handle empty or very small ranges
        if (total_elements == 0) {
            return std::vector<math::Range>(B, {range.start, range.start});
        }

        if (B <= 1 && !force_B.has_value()) {
            return {range};
        }

        const size_t vectors_per_block = total_packed / B;
        const size_t remainder_vectors = total_packed % B;

        std::vector<math::Range> blocks(B);
        size_t current = range.start;

        for (size_t i = 0; i < B; ++i) {
            const size_t v_count = vectors_per_block + (i < remainder_vectors ? 1 : 0);
            size_t size = v_count * v_size;

            // Attach scalar tail to the very last block
            if (i == B - 1) {
                size += total_elements % v_size;
            }

            blocks[i] = {current, current + size};
            current += size;
        }

        return blocks;
    }


    // partition a stl container into equally sized work spans
    template <std::ranges::contiguous_range R>
    [[nodiscard]] auto make_linear_schedule(
        R&& container,
        const BlockConfig& config = {})
    {
        using T = std::ranges::range_value_t<R>;

        // Generate aligned integer ranges
        const math::Range full_range{0, std::ranges::size(container)};
        const auto block_ranges = make_linear_schedule(full_range, config);

        // Map integer ranges to std::span slices
        std::vector<std::span<T>> spans;
        spans.reserve(block_ranges.size());

        auto* data_ptr = std::ranges::data(container);
        for (const auto& r : block_ranges) {
            spans.emplace_back(data_ptr + r.start, r.size());
        }

        return spans;
    }



    // Schedule pairwise self-interactions into independent execution phases
    [[nodiscard]] inline SymmetricPhaseSchedule make_symmetric_schedule(
        const math::Range& range,
        const BlockConfig& config)
    {
        SymmetricPhaseSchedule schedule;
        if (range.empty()) return schedule;

        // 1. Calculate block count (B)
        const size_t max_aligned_units = range.size() / config.alignment;
        size_t B = std::max<size_t>(1, config.target_chunk_size);

        if (B > max_aligned_units && max_aligned_units > 0) {
            B = max_aligned_units;
        }

        // 2. Partition range
        std::vector<math::Range> blocks = make_linear_schedule(range, config);
        B = blocks.size();

        // 3. Force even B for tournament math by adding a dummy bye block
        if (B % 2 != 0) {
            blocks.push_back({range.stop, range.stop});
            B++;
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
    [[nodiscard]] inline BipartitePhaseSchedule make_bipartite_schedule(
        const math::Range& range1,
        const math::Range& range2,
        const BlockConfig& config)
    {
        BipartitePhaseSchedule schedule;
        if (range1.empty() || range2.empty()) return schedule;

        // 1. Determine common block count B for a square task matrix
        size_t B = std::max<size_t>(1, config.target_chunk_size);

        const size_t max_aligned1 = range1.size() / config.alignment;
        const size_t max_aligned2 = range2.size() / config.alignment;
        const size_t max_aligned = std::max(max_aligned1, max_aligned2);

        if (B > max_aligned && max_aligned > 0) {
            B = max_aligned;
        }

        // 2. Partition both ranges into exactly B blocks
        std::vector<math::Range> blocks1 = make_linear_schedule(range1, config, B);
        std::vector<math::Range> blocks2 = make_linear_schedule(range2, config, B);

        AP_ASSERT(blocks1.size() == blocks2.size(), "");

        // 3. Cyclic Diagonal Scheduling (Phase p: Pair block i with block (i+p)%B)
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