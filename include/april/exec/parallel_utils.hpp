#pragma once

#include <vector>
#include <algorithm>
#include <span>
#include <ranges>

#include "april/exec/info.hpp"
#include "april/base/types.hpp"
#include "april/math/range.hpp"

namespace april::exec {

    struct BlockConfig {
        size_t num_threads = N_CPU_THREADS;
        size_t target_chunk_size = 256;
        size_t oversubscription = 4;
        size_t alignment = packed::size();
    };

    // partition a range into a bunch of equal sized chunks
    [[nodiscard]] inline std::vector<math::Range> partition_work(
        const math::Range& range,
        const BlockConfig& config = {})
    {
        const size_t total_elements = range.size();
        if (total_elements == 0) return {};

        const size_t v_size = config.alignment;
        const size_t total_packed = total_elements / v_size;

        // Determine target block count
        size_t B = std::max<size_t>(1, config.num_threads * config.oversubscription);

        // Clamp chunks to avoid starving threads of full SIMD vectors
        if (total_packed < B) {
            B = std::max<size_t>(1, total_packed);
        }

        // Fast path for sequential/small execution
        if (B <= 1) {
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


    // partition a stl container into a bunch of equal sized chunks
    template <std::ranges::contiguous_range R>
    [[nodiscard]] auto partition_work(
        R&& container,
        const BlockConfig& config = {})
    {
        using T = std::ranges::range_value_t<R>;

        // Generate aligned integer ranges
        const math::Range full_range{0, std::ranges::size(container)};
        const auto block_ranges = partition_work(full_range, config);

        // Map integer ranges to std::span slices
        std::vector<std::span<T>> spans;
        spans.reserve(block_ranges.size());

        auto* data_ptr = std::ranges::data(container);
        for (const auto& r : block_ranges) {
            spans.emplace_back(data_ptr + r.start, r.size());
        }

        return spans;
    }

} // namespace april::exec