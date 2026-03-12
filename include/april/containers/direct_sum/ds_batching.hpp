#pragma once

#include <barrier>
#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"
#include "april/exec/info.hpp"


#include "april/exec/executors/executor_traits.hpp"


namespace april::container::internal {
    template <typename Container, typename AsymmetricBatch, typename SymmetricBatch, exec::IsExecutor Executor>
    struct SymmetricParallelBatch : batching::BatchBase<exec::ParallelTrait::IntraBatch,
                                                        AsymmetricBatch::vector_trait & SymmetricBatch::vector_trait> {
        explicit SymmetricParallelBatch(Container& container, const Executor& executor)
            : container(container), executor(executor) {}

        template <ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            executor.execute(diagonal_phase.size(), [&](size_t i) {
                diagonal_phase[i].template for_each_pair<P, E>(f);
            });

            for (auto& phase : off_diagonal_phases) {
                executor.execute(phase.size(), [&](size_t i) {
                    phase[i].template for_each_pair<P, E>(f);
                });
            }
        }


        void set_range(const math::Range& range, const size_t oversubscription = 4) {
            const size_t n_threads = exec::N_CPU_THREADS;
            constexpr size_t target_chunk_size = 256;
            size_t B = std::max<size_t>(1, range.size() / target_chunk_size);

            // Make sure B is even for the tournament math
            if (B % 2 != 0) B++;

            // Also ensure B is large enough to keep threads busy (oversubscription check)
            const size_t min_B = n_threads * oversubscription;
            if (B < min_B) B = (min_B % 2 == 0) ? min_B : min_B + 1;

            // Partition the range into B blocks with strict SIMD alignment
            constexpr size_t v_size = packed::size();
            const size_t total_vectors = range.size() / v_size;
            const size_t vectors_per_block = total_vectors / B;
            const size_t remainder_vectors = total_vectors % B;

            std::vector<math::Range> blocks(B);
            size_t current_idx = range.start;

            for (size_t i = 0; i < B; ++i) {
                const size_t v_count = vectors_per_block + (i < remainder_vectors ? 1 : 0);
                size_t size = v_count * v_size;

                // The final block picks up any trailing scalar elements
                if (i == B - 1) {
                    size += range.size() % v_size;
                }

                blocks[i] = {current_idx, current_idx + size};
                current_idx += size;
            }

            // Phase 0: create symmetric batches for the diagonals
            for (size_t i = 0; i < B; ++i) {
                if (blocks[i].start == blocks[i].stop) continue;

                SymmetricBatch sym(this->container);
                sym.types = this->types;
                sym.range = blocks[i];
                diagonal_phase.push_back(sym);
            }

            // Phases 1 to B-1: create asymmetric batches for off diagonals
            off_diagonal_phases.resize(B - 1);
            std::vector<size_t> circle(B);
            for (size_t i = 0; i < circle.size(); ++i) circle[i] = i;

            for (size_t phase = 0; phase < B - 1; ++phase) {
                for (size_t i = 0; i < B / 2; ++i) {
                    // block b1 interacts with block b2
                    size_t b1 = circle[i];
                    size_t b2 = circle[B - 1 - i];

                    const size_t idx1 = std::min(b1, b2);
                    const size_t idx2 = std::max(b1, b2);

                    // if empty skip
                    if (blocks[idx1].start == blocks[idx1].stop ||
                        blocks[idx2].start == blocks[idx2].stop)
                        continue;

                    AsymmetricBatch asym(container);
                    asym.types = this->types;
                    asym.range1 = blocks[idx1];
                    asym.range2 = blocks[idx2];
                    off_diagonal_phases[phase].push_back(asym);
                }

                // rotate circle array for next phase: keep index 0 fixed, rotate 1..B-1 right
                const size_t last = circle.back();
                for (size_t i = B - 1; i > 1; --i) {
                    circle[i] = circle[i - 1];
                }
                circle[1] = last;
            }
        }

    private:
        Container& container;
        const Executor& executor;
        std::vector<SymmetricBatch> diagonal_phase;
        std::vector<std::vector<AsymmetricBatch>> off_diagonal_phases;
    };


    template <typename Container, typename AsymmetricBatch, exec::IsExecutor Executor>
    struct
        AsymmetricParallelBatch : batching::BatchBase<exec::ParallelTrait::IntraBatch, AsymmetricBatch::vector_trait> {
        explicit AsymmetricParallelBatch(Container& container, const Executor& executor)
            : container(container), executor(executor) {}

        template <ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            // The asymmetric case is strictly bipartite, so all phases are structurally
            // identical. There is no separate diagonal phase to handle.
            for (const auto& phase : phases) {
                executor.execute(phase.size(), [&](size_t i) {
                    phase[i].template for_each_pair<P, E>(f);
                });
            }
        }

        void set_range(const math::Range& range1, const math::Range& range2, const size_t oversubscription = 4) {
            const size_t n_threads = exec::N_CPU_THREADS;
            constexpr size_t target_chunk_size = 256;

            // Use the larger range to determine block count
            const size_t max_range = std::max(range1.size(), range2.size());
            size_t B = std::max<size_t>(1, max_range / target_chunk_size);

            const size_t min_B = n_threads * oversubscription;
            if (B < min_B) B = min_B;

            // Partition BOTH ranges using SIMD-aligned boundaries
            std::vector<math::Range> blocks1 = partition_range(range1, B);
            std::vector<math::Range> blocks2 = partition_range(range2, B);

            // 3. Schedule using Cyclic Diagonals
            phases.resize(B);
            for (size_t p = 0; p < B; ++p) {
                for (size_t i = 0; i < B; ++i) {
                    // Modulo shift ensures no column overlap within the phase
                    const size_t j = (i + p) % B;

                    // Skip if either block ended up empty (e.g., if B > total particles)
                    if (blocks1[i].start == blocks1[i].stop ||
                        blocks2[j].start == blocks2[j].stop) {
                        continue;
                    }

                    AsymmetricBatch asym(container);
                    asym.types = this->types;
                    asym.range1 = blocks1[i];
                    asym.range2 = blocks2[j];
                    phases[p].push_back(asym);
                }
            }
        }

        std::vector<std::vector<AsymmetricBatch>> phases;

    private:
        Container& container;
        const Executor& executor;

        // Helper to divide a range into exactly B chunks, distributing the remainder
        static std::vector<math::Range> partition_range(const math::Range& r, const size_t num_blocks) {
            constexpr size_t v_size = packed::size();
            const size_t total_vectors = r.size() / v_size;
            const size_t vectors_per_block = total_vectors / num_blocks;
            const size_t remainder_vectors = total_vectors % num_blocks;

            std::vector<math::Range> blocks(num_blocks);
            size_t current = r.start;

            for (size_t i = 0; i < num_blocks; ++i) {
                const size_t v_count = vectors_per_block + (i < remainder_vectors ? 1 : 0);
                size_t size = v_count * v_size;

                // The final block picks up any trailing scalar elements
                if (i == num_blocks - 1) {
                    size += r.size() % v_size;
                }

                blocks[i] = {current, current + size};
                current += size;
            }
            return blocks;
        }
    };


    template <typename Container, typename AsymmetricBatch, typename SymmetricBatch, typename ChunkPtr, exec::IsExecutor
              Executor>
    struct SymmetricParallelChunkedBatch : batching::BatchBase<exec::ParallelTrait::IntraBatch,
                                                               AsymmetricBatch::vector_trait &
                                                               SymmetricBatch::vector_trait> {
        // Injected ChunkPtr* ensures we don't need to access protected members of Container directly
        explicit SymmetricParallelChunkedBatch(Container& container, ChunkPtr* chunks, const Executor& executor)
            : container(container), chunks(chunks), executor(executor) {}

        template <ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            executor.execute(diagonal_phase.size(), [&](size_t i) {
                diagonal_phase[i].template for_each_pair<P, E>(f);
            });

            for (const auto& phase : off_diagonal_phases) {
                executor.execute(phase.size(), [&](size_t i) {
                    phase[i].template for_each_pair<P, E>(f);
                });
            }
        }

        void set_range(const math::Range& range_chunks, const size_t tail, const size_t oversubscription = 4) {
            const size_t n_threads = executor.num_threads();

            size_t B = std::max<size_t>(1, n_threads * oversubscription);
            if (B % 2 != 0) B++;

            if (B > range_chunks.size() && range_chunks.size() > 0) {
                B = range_chunks.size();
                if (B % 2 != 0) B--;
                if (B == 0) B = 2;
            }

            // 1. Partition the chunk-indices into B blocks
            std::vector<math::Range> blocks(B);
            const size_t chunks_per_block = range_chunks.size() / B;
            const size_t remainder = range_chunks.size() % B;
            size_t current_idx = range_chunks.start;

            for (size_t i = 0; i < B; ++i) {
                const size_t size = chunks_per_block + (i < remainder ? 1 : 0);
                blocks[i] = {current_idx, current_idx + size};
                current_idx += size;
            }

            // 2. Phase 0: Diagonals
            for (size_t i = 0; i < B; ++i) {
                if (blocks[i].empty()) continue;

                // Pass the injected chunks pointer to the underlying batch
                SymmetricBatch sym(this->container, this->chunks);
                sym.types = this->types;
                sym.range_chunks = blocks[i];
                sym.range_tail = (blocks[i].stop == range_chunks.stop) ? tail : 0;
                diagonal_phase.push_back(std::move(sym));
            }

            // 3. Phases 1 to B-1: Off-diagonals (Tournament)
            off_diagonal_phases.resize(B - 1);
            std::vector<size_t> circle(B);
            for (size_t i = 0; i < circle.size(); i ++) circle[i] = i;

            for (size_t phase = 0; phase < B - 1; ++phase) {
                for (size_t i = 0; i < B / 2; ++i) {
                    size_t b1 = circle[i];
                    size_t b2 = circle[B - 1 - i];

                    const size_t idx1 = std::min(b1, b2);
                    const size_t idx2 = std::max(b1, b2);

                    if (blocks[idx1].empty() || blocks[idx2].empty()) continue;

                    // Pass the injected chunks pointer to the underlying batch
                    AsymmetricBatch asym(this->container, this->chunks);
                    asym.types = this->types;
                    asym.range1_chunks = blocks[idx1];
                    asym.range2_chunks = blocks[idx2];

                    asym.range1_tail = (blocks[idx1].stop == range_chunks.stop) ? tail : 0;
                    asym.range2_tail = (blocks[idx2].stop == range_chunks.stop) ? tail : 0;

                    off_diagonal_phases[phase].push_back(std::move(asym));
                }

                size_t last = circle.back();
                for (size_t i = B - 1; i > 1; --i) circle[i] = circle[i - 1];
                circle[1] = last;
            }
        }

    private:
        Container& container;
        ChunkPtr* chunks; // Stored pointer
        const Executor& executor;
        std::vector<SymmetricBatch> diagonal_phase;
        std::vector<std::vector<AsymmetricBatch>> off_diagonal_phases;
    };


    template <typename Container, typename AsymmetricBatch, typename ChunkPtr, exec::IsExecutor Executor>
    struct AsymmetricParallelChunkedBatch : batching::BatchBase<exec::ParallelTrait::IntraBatch,
                                                                AsymmetricBatch::vector_trait> {
        explicit AsymmetricParallelChunkedBatch(Container& container, ChunkPtr* chunks, const Executor& executor)
            : container(container), chunks(chunks), executor(executor) {}

        template <ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            for (const auto& phase : phases) {
                executor.execute(phase.size(), [&](size_t i) {
                    phase[i].template for_each_pair<P, E>(f);
                });
            }
        }

        void set_range(const math::Range& range1_chunks, size_t tail1,
                       const math::Range& range2_chunks, size_t tail2,
                       const size_t oversubscription = 4) {
            const size_t n_threads = executor.num_threads();
            const size_t max_chunks = std::max(range1_chunks.size(), range2_chunks.size());

            // 1. Determine common block count B for a square task matrix
            size_t B = std::max<size_t>(1, n_threads * oversubscription);

            // Clamp B if it significantly exceeds the available chunks to avoid excessive empty tasks
            if (B > max_chunks && max_chunks > 0) B = max_chunks;

            // 2. Partition both ranges into exactly B blocks
            auto blocks1 = partition_chunks(range1_chunks, B);
            auto blocks2 = partition_chunks(range2_chunks, B);

            // 3. Schedule using Cyclic Diagonals (Phase p: Pair block i with block (i+p)%B)
            phases.resize(B);
            for (size_t p = 0; p < B; ++p) {
                for (size_t i = 0; i < B; ++i) {
                    size_t j = (i + p) % B;

                    if (blocks1[i].empty() || blocks2[j].empty()) continue;

                    AsymmetricBatch asym(container, chunks);
                    asym.types = this->types;

                    asym.range1_chunks = blocks1[i];
                    asym.range2_chunks = blocks2[j];

                    // Apply tails only if the sub-block reaches the end of the original range
                    asym.range1_tail = (blocks1[i].stop == range1_chunks.stop) ? tail1 : 0;
                    asym.range2_tail = (blocks2[j].stop == range2_chunks.stop) ? tail2 : 0;

                    phases[p].push_back(std::move(asym));
                }
            }
        }

    private:
        Container& container;
        ChunkPtr* chunks;
        const Executor& executor;
        std::vector<std::vector<AsymmetricBatch>> phases;

        // Internal helper to force a range into exactly B chunks
        static std::vector<math::Range> partition_chunks(const math::Range& r, size_t B) {
            std::vector<math::Range> res(B);
            if (r.empty()) return res;

            const size_t size = r.size() / B;
            const size_t rem = r.size() % B;
            size_t curr = r.start;

            for (size_t i = 0; i < B; ++i) {
                const size_t count = size + (i < rem ? 1 : 0);
                res[i] = {curr, curr + count};
                curr += count;
            }
            return res;
        }
    };
}
