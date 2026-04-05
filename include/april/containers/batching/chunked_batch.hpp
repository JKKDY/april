/**
* @file chunked_batch.hpp
 * @brief Core interaction code for AoSoA memory layouts. Orchestrates SIMD and scalar computation.
 * * This file implements the Rotation Sweep algorithm for full N x N pair interactions directly
 * within registers. It uses a three-tier partitioning strategy on SIMD execution:
 * 1. Body: Full chunks processed via unmasked SIMD.
 * 2. Full Tail: SIMD-aligned blocks within partial chunks (also unmasked).
 * 3. Partial Tail: Scalar remainders handled via masking and broadcasting.
 */

#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"


namespace april::container::batching {

    namespace internal {

        /**
         * 360-degree SIMD rotation sweep for distinct blocks (N^2 within the blocks).
         * Pairs every lane in b1 with every lane in acc2 (accumulator) via cyclic rotation.
         * Reciprocal forces are accumulated into acc2 memory at the end.
         */
        template <typename Buffer1, typename PackedAccessor, typename Kernel>
        APRIL_FORCE_INLINE inline void interact_block_vs_block(Buffer1& b1, PackedAccessor&& acc2, Kernel& f) {
            auto b2 = acc2.load_buffer();

            // interact all pairs {p1, p2} with p1 in b1, p2 in b2
            APRIL_UNROLL_LOOP_N(packed::size())
            for (size_t k = 0; k < packed::size(); k++) {
                f(b1.to_view(), b2.to_view());
                b2.rotate_right();
            }

            b2.update_into(acc2); // Flush accumulated forces back to memory
        }


        /**
         * Compile-time resolver for container data access.
         * Simplifies batch loops by wrapping the container's direct memory accessors
         * with the specific Read/Write masks required by the kernel.
         */
        template<typename Container, typename Kernel>
        struct BatchContext {
            using K = std::remove_cvref_t<Kernel>;
            Container& container;
            K & kernel;

            BatchContext(Container& c, K& k) : container(c), kernel(k) {}

            APRIL_FORCE_INLINE auto scalar(size_t chunk, size_t i) const {
                return container.template at<K::Read, K::Write>(chunk, i);
            }
            APRIL_FORCE_INLINE auto packed(size_t chunk, size_t i) const {
                return container.template at_packed<K::Read, K::Write>(chunk, i);
            }

            APRIL_FORCE_INLINE void prefetch(size_t, size_t = 0) const {
                // container.template prefetch_particle<K::Read | K::Write>(chunk, i);
            }

            APRIL_FORCE_INLINE void prefetch_nta(size_t , size_t = 0) const {
                // container.template prefetch_particle_nta<K::Read | K::Write>(chunk, i);
            }
        };
        // Deduction guide for convenience
        template <typename Container, typename Kernel>
        BatchContext(Container&, Kernel&) -> BatchContext<Container, Kernel>;


        /**
         * Base class for chunk-based (AoSoA) batch iterators.
         * Handles shared constants and the static lane-index array used for SIMD masking
         * in partial chunks.
         */
        template <typename Container, typename ChunkPtr>
        struct ChunkedBatchBase : BatchBase<exec::VectorTrait::ScalarPath | exec::VectorTrait::VectorPath> {
            explicit ChunkedBatchBase(Container& container, ChunkPtr* chunks)
                : container(container), chunks(chunks) {
                for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
            }

            /**
             * Dispatch for pair-wise iteration. Applies a kernel f to each pair
             * Branches to the vectorized path if the policy allows and the kernel supports it.
             */
            template <exec::ExecutionMode E, exec::IsKernel Func>
            void for_each_pair(this const auto & self, Func&& f) {
                if (self.empty()) return;

                if constexpr (static_cast<bool>(E & exec::ExecutionMode::Vector))
                    self.for_each_pair_packed(std::forward<Func>(f));
                else
                    self.for_each_pair_scalar(std::forward<Func>(f));
            }

        protected:
            Container& container;
            ChunkPtr* APRIL_RESTRICT const chunks;
            static constexpr size_t chunk_size = Container::chunk_size; // power of 2
            static constexpr size_t packed_size = packed::size(); // power of 2
            static constexpr size_t iter_chunks = chunk_size / packed_size; // how many packed units in a chunk
            alignas(64) packed::value_type idx_arr[packed_size]{};
        };
    }


    //-----------------
    // ASYMMETRIC BATCH
    //-----------------
    /**
     * Batch for interacting two bipartite particle sets of particles
     */
    template <typename Container, typename ChunkPtr>
    struct AsymmetricChunkedBatch : internal::ChunkedBatchBase<Container, ChunkPtr> {
        using Base = internal::ChunkedBatchBase<Container, ChunkPtr>;
        using Base::Base, Base::container, Base::chunks, Base::chunk_size, Base::packed_size, Base::iter_chunks, Base::idx_arr;
        friend Base;

        [[nodiscard]] bool empty() const noexcept{
            return range1_chunks.empty() || range2_chunks.empty();
        }

        // range1_chunks/range2_chunks represent chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
        math::Range range1_chunks; // Chunk Indices [start, end)
        math::Range range2_chunks;
        size_t range1_tail{}; // Number of valid items in the last chunk of range1 (0 = Full)
        size_t range2_tail{};

    private:

        //------------
        // SCALAR PATH
        //------------
        /**
         * Scalar path for asymmetric interactions.
         * Simple nested loops covering all combinations, including tail remainders.
         */
        template <exec::IsKernel Kernel>
        void for_each_pair_scalar(Kernel && f) const {
            internal::BatchContext ctx(container, f);

            // peel of last chunk (i.e. the tail)
            const size_t c1_body_end = range1_chunks.stop - 1;
            const size_t c2_body_end = range2_chunks.stop - 1;

            // if tail is 0, it means the chunk is actually full (stride)
            const size_t limit1_tail = (range1_tail == 0) ? chunk_size : range1_tail;
            const size_t limit2_tail = (range2_tail == 0) ? chunk_size : range2_tail;

            // body vs body
            for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
                ctx.prefetch( c1 + 1);

                for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
                    ctx.prefetch_nta( c2 + 1);

                    for (size_t i = 0; i < chunk_size; ++i) {
                        auto p1 = ctx.scalar(c1, i);
                        for (size_t j = 0; j < chunk_size; ++j) {
                            auto p2 = ctx.scalar(c2, j);
                            f(p1, p2);
                        }
                    }
                }
            }

            // body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
            for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
                ctx.prefetch( c1 + 1);
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = ctx.scalar(c1, i);
                    for (size_t j = 0; j < limit2_tail; ++j) {
                        auto p2 = ctx.scalar(c2_body_end, j);
                        f(p1, p2);
                    }
                }
            }

            // body 2 vs tail 1 (iterate full chunks of R2 against the single partial chunk of R1
            for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
                ctx.prefetch( c2 + 1);
                for (size_t i = 0; i < limit1_tail; ++i) {
                    auto p1 = ctx.scalar(c1_body_end, i);
                    for (size_t j = 0; j < chunk_size; ++j) {
                        auto p2 = ctx.scalar(c2, j);
                        f(p1, p2);
                    }
                }
            }

            // tail 1 vs tail 2 (Interaction between the two last chunks)
            for (size_t i = 0; i < limit1_tail; ++i) {
                auto p1 = ctx.scalar(c1_body_end, i);
                for (size_t j = 0; j < limit2_tail; ++j) {
                    auto p2 = ctx.scalar(c2_body_end, j);
                    f(p1, p2);
                }
            }
        }


        //----------
        // SIMD PATH
        //----------
        /**
         * Main vectorized path for asymmetric interactions.
         */
        template <exec::IsKernel Kernel>
        void for_each_pair_packed(Kernel&& f) const {
            // handle full chunks masquerading as 0-length tails
            const size_t tail1_end = (range1_tail == 0) ? chunk_size : range1_tail;
            const size_t tail2_end = (range2_tail == 0) ? chunk_size : range2_tail;

            // calculate exact index limits of full tail (largest index that is a multiple of packed_size)
            const size_t full_tail_end1 = (tail1_end / packed_size) * packed_size;
            const size_t full_tail_end2 = (tail2_end / packed_size) * packed_size;

            // we split set1 (range1_chunks + range1_tail) and set2 (range2_chunks + range2_tail) into 3 parts each:
            // Chunk indices excluding the final (potentially partial) tail chunk
            const math::Range full_chunks1 = {range1_chunks.start, range1_chunks.stop - 1};
            const math::Range full_chunks2 = {range2_chunks.start, range2_chunks.stop - 1};

            // SIMD-aligned offsets within the tail chunk (increments by packed_size)
            const math::Range full_tail1 = {0, full_tail_end1, packed_size};
            const math::Range full_tail2 = {0, full_tail_end2, packed_size};

            // Scalar remainder indices at the end of the tail chunk (single particle steps)
            const math::Range partial_tail1 = {full_tail_end1, tail1_end};
            const math::Range partial_tail2 = {full_tail_end2, tail2_end};

            // we interact them as follows (set1 = [C1, F1, P1], set2 = [C2, F2, P2]):
            // 1. [C1] vs [a) C2 + b) F2]  --> Body1 SIMD vs Body2 + Tail2 SIMD
            // 2. [F1] vs [a) C2 + b) F2]  --> Tail1 SIMD vs Body2 + Tail2 SIMD
            interact_body_vs_body(full_chunks1, full_chunks2, full_tail1, full_tail2, f);

            // 3. [P2] vs [a) C1 + b) F1]  --> Tail2 Scalar vs Body1 + Tail1 SIMD
            // 4. [P1] vs [a) C2 + b) F2]  --> Tail1 Scalar vs Body2 + Tail2 SIMD
            interact_tails_vs_body(full_chunks1, full_chunks2, full_tail1, full_tail2, partial_tail1, partial_tail2, f);

            // 5. [P1] vs [P2]  --> Scalar vs Scalar Remainder
            interact_tails_vs_tails(full_chunks1, full_chunks2, partial_tail1, partial_tail2, f);
        }


        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_body_vs_body(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& full_tail1, const math::Range& full_tail2,
            Kernel && f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 1. full SIMD blocks in range 1 (Chunks) vs full SIMD blocks in range 2 (Chunks + Tail)
            for (size_t c1 : full_chunks1) {
                ctx.prefetch( c1 + 1);
                APRIL_UNROLL_LOOP_N(iter_chunks)
                for (size_t i = 0; i < chunk_size; i += packed_size) {
                    auto packed1 = ctx.packed(c1, i);
                    auto buffer1 = packed1.load_buffer(); // load simd block in range 1

                    // a. sweep buffer1 across all full SIMD blocks in Range 2 body chunks [C2]
                    for (size_t c2 : full_chunks2) {
                        ctx.prefetch_nta( c2 + 1);
                        APRIL_UNROLL_LOOP_N(iter_chunks)
                        for (size_t j = 0; j < chunk_size; j += packed_size) {
                            interact_block_vs_block(buffer1, ctx.packed(c2, j), f);
                        }
                    }

                    // b. sweep buffer1 across SIMD-aligned blocks within the Range 2 tail chunk [F2]
                    for (size_t t2 : full_tail2) {
                        interact_block_vs_block(buffer1, ctx.packed(full_chunks2.stop, t2), f);
                    }

                    buffer1.update_into(packed1); // final register flush for buffer1
                }
            }

            // 2. SIMD-aligned blocks in the Range 1 tail chunk vs Full Range 2 (Body + Tail)
            for (size_t t1 : full_tail1) {
                auto packed1 = ctx.packed(full_chunks1.stop, t1);
                auto buffer1 = packed1.load_buffer();

                // a. Interaction with all full chunks in the Range 2 body [C2]
                for (size_t c2 : full_chunks2) {
                    ctx.prefetch_nta( c2 + 1);
                    APRIL_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size) {
                        interact_block_vs_block(buffer1, ctx.packed(c2, j), f);
                    }
                }

                // b. Interaction with SIMD-aligned blocks within the Range 2 tail chunk [F2]
                for (size_t t2 : full_tail2) {
                    interact_block_vs_block(buffer1, ctx.packed(full_chunks2.stop, t2), f);
                }
                buffer1.update_into(packed1);
            }
        }


         // Range 1 is the SIMD Block, Range 2 is the Broadcast Scalar
        template <typename PackedAccessor, typename BufferScalar, typename Kernel>
        APRIL_FORCE_INLINE void interact_block1_vs_scalar2(PackedAccessor&& p1_block, BufferScalar& p2_scalar, Kernel& f) const {
            auto b_block = p1_block.load_buffer();
            f(b_block.to_view(), p2_scalar.to_view()); // P1 first, P2 second
            b_block.update_into(p1_block);
        }


        // Range 1 is the Broadcast Scalar, Range 2 is the SIMD Block
        template <typename BufferScalar, typename PackedAccessor, typename Kernel>
        APRIL_FORCE_INLINE void interact_scalar1_vs_block2(BufferScalar& p1_scalar, PackedAccessor&& p2_block, Kernel& f) const {
            auto b_block = p2_block.load_buffer();
            f(p1_scalar.to_view(), b_block.to_view()); // P1 first, P2 second
            b_block.update_into(p2_block);
        }


        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_tails_vs_body(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& full_tail1,   const math::Range& full_tail2,
            const math::Range& partial_tail1,const math::Range& partial_tail2,
            Kernel&& f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 3. Range 2 Partial Tail [P2] vs Range 1 SIMD blocks [C1 + F1]
            for (size_t i : partial_tail2) {
                auto p2 = ctx.scalar(full_chunks2.stop, i);
                auto buffer2 = p2.broadcast();

                // a. Interaction with full chunks in Range 1 body
                for (size_t c1 : full_chunks1) {
                    ctx.prefetch_nta(c1 + 1);
                    APRIL_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size)
                        interact_block1_vs_scalar2(ctx.packed(c1, j), buffer2, f);
                }

                // b. Interaction with SIMD-aligned blocks in Range 1 tail chunk
                for (size_t t1 : full_tail1)
                    interact_block1_vs_scalar2(ctx.packed(full_chunks1.stop, t1), buffer2, f);

                buffer2.reduce_into(p2);
            }

            // Partial Tail 1 vs Full Range 2 (chunks + full tail)
            for (size_t i : partial_tail1) {
                auto p1 = ctx.scalar(full_chunks1.stop, i);
                auto buffer1 = p1.broadcast();

                // a. Interaction with full chunks in Range 2 body
                for (size_t c2 : full_chunks2) {
                    ctx.prefetch_nta( c2 + 1);
                    APRIL_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size)
                        interact_scalar1_vs_block2(buffer1, ctx.packed(c2, j), f);
                }

                // b. Interaction with SIMD-aligned blocks in Range 2 tail chunk
                for (size_t t2 : full_tail2)
                    interact_scalar1_vs_block2(buffer1, ctx.packed(full_chunks2.stop, t2), f);

                buffer1.reduce_into(p1);
            }
        }

        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_tails_vs_tails(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& partial_tail1, const math::Range& partial_tail2,
            Kernel && f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 5. partial tail vs partial tail
            if (partial_tail1.start != partial_tail1.stop && partial_tail2.start != partial_tail2.stop) {
                // initialize mask for the valid lanes in partial_tail1
                const auto lane_indices = packed::load_aligned(idx_arr);
                const auto valid_lanes = static_cast<double>(partial_tail1.stop - partial_tail1.start);
                const auto mask = lane_indices < valid_lanes;

                // load the single SIMD block containing partial_tail1
                auto packed1 = ctx.packed(full_chunks1.stop, partial_tail1.start);
                auto buffer1 = packed1.load_buffer();

                // loop over the individual particles in partial_tail2
                for (size_t i = partial_tail2.start; i < partial_tail2.stop; ++i) {
                    auto p2 = ctx.scalar(full_chunks2.stop, i);
                    auto buffer2 = p2.broadcast();

                    auto view1 = buffer1.to_view();
                    auto view2 = buffer2.to_view();
                    f(view1, view2);

                    buffer2.reduce_into(p2, mask);
                }

                // write back valid lanes to memory for packed1
                buffer1.update_into(packed1, mask);
            }
        }
    };





    //----------------
    // SYMMETRIC BATCH
    //----------------
    /**
     * Batch for self-interaction of a set of particles.
     * Calculates the upper triangle (j > i) to avoid double-counting.
     * Uses a specialized 180-degree rotation trick for Register-level symmetry.
     */
    template <typename Container, typename ChunkPtr>
    struct SymmetricChunkedBatch : internal::ChunkedBatchBase<Container, ChunkPtr> {
        using Base = internal::ChunkedBatchBase<Container, ChunkPtr>;
        using Base::Base, Base::container, Base::chunks, Base::chunk_size, Base::packed_size, Base::iter_chunks, Base::idx_arr;
        friend Base;

        [[nodiscard]] bool empty() const noexcept{
            return range_chunks.start == range_chunks.stop;
        }

        // Range represents chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
        math::Range range_chunks; // Chunk Indices [start, end)
        size_t range_tail{}; // Number of valid items in the last chunk of range1 (0 = Full)
    private:

        //----------------
        // SCALAR ITERATOR
        //----------------
        template <typename Kernel>
        void for_each_pair_scalar(Kernel && f) const {
            internal::BatchContext ctx(container, f);

            const size_t c_body_end = range_chunks.stop - 1;
            // Note: Use chunk_size here for full chunks, not packed_size
            const size_t limit_tail = (range_tail == 0) ? chunk_size : range_tail;

            for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
                ctx.prefetch( c1 + 1);

                // chunk self interaction
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = ctx.scalar(c1, i);
                    for (size_t j = i + 1; j < chunk_size; ++j) {
                        auto p2 = ctx.scalar(c1, j);
                        f(p1, p2);
                    }
                }

                // interaction with all other chunks
                for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
                    ctx.prefetch_nta( c2 + 1);
                    for (size_t i = 0; i < chunk_size; ++i) {
                        auto p1 = ctx.scalar(c1, i);
                        for (size_t j = 0; j < chunk_size; ++j) {
                            auto p2 = ctx.scalar(c2, j);
                            f(p1, p2);
                        }
                    }
                }
            }

            // body vs tail (every body chunk with tail chunk)
            for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
                ctx.prefetch( c1 + 1);
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = ctx.scalar(c1, i);
                    for (size_t j = 0; j < limit_tail; ++j) {
                        auto p2 = ctx.scalar(c_body_end, j);
                        f(p1, p2);
                    }
                }
            }

            // tail (interact tail chunk with itself)
            for (size_t i = 0; i < limit_tail; ++i) {
                auto p1 = ctx.scalar(c_body_end, i);
                for (size_t j = i + 1; j < limit_tail; ++j) {
                    auto p2 = ctx.scalar(c_body_end, j);
                    f(p1, p2);
                }
            }
        }


        //----------------
        // PACKED ITERATOR
        //----------------
        /**
         * Symmetric self-interaction sweep within a SIMD block.
         * Performs N/2 rotations. For N=8, rotates 3 times to cover 24 pairs,
         * then realigns and performs a final 180-degree flip to capture the last 4 pairs.
         * This ensures every unique pair (i, j) is calculated exactly once.
         */
        template <typename BufferT, typename PackedAccessor, typename Kernel>
        APRIL_FORCE_INLINE void interact_symmetric_self(BufferT& buffer1, const PackedAccessor& packed1, Kernel& f) const {
            auto buffer2 = packed1.load_buffer();

            APRIL_UNROLL_LOOP()
            for (size_t k = 0; k < packed_size / 2 - 1; k++) {
                buffer2.rotate_right();
                auto view1 = buffer1.to_view();
                auto view2 = buffer2.to_view();
                f(view1, view2);
            }

            buffer2.template rotate_left<packed_size / 2 - 1>(); // Realign buffer2 back to its origin lanes
            buffer1.accumulate(buffer2); // merge deltas from buffer 2 into buffer 1
            buffer2.template rotate_right<packed_size / 2>(); // rotate 180 for last half step

            auto view1 = buffer1.to_view();
            auto view2 = buffer2.to_view();
            f(view1, view2);
        }

        // helper for packed vs broadcast scalar block (no rotations needed)
        template <typename BufferT, typename ScalarAccessor, typename Kernel>
        APRIL_FORCE_INLINE void interact_block_vs_scalar(BufferT& buffer_block, const ScalarAccessor& p2, Kernel& f) const {
            auto buffer_scalar = p2.broadcast();
            f(buffer_block.to_view(), buffer_scalar.to_view());
            buffer_scalar.reduce_into(p2);
        }

        /**
         * Main vectorized path for symmetric interactions.
         * Covers:
         * 1. Block Self-interaction (interact_symmetric_self).
         * 2. Intra-chunk/Inter-chunk block pairs (interact_block_vs_block).
         * 3. Remainder tail masking via upper-triangle indices.
         */
        void for_each_pair_packed(auto&& f) const {
            const size_t tail_len = (range_tail == 0) ? chunk_size : range_tail; // length of entire tail
            const size_t full_tail_end = (tail_len / packed_size) * packed_size;
            // largest index that is a multiple of packed_size

            // Partitioning the chunked container into 3 distinct sets:
            // 1. C: Full chunks (excluding the final tail chunk)
            // 2. F: Full SIMD blocks within the tail chunk (SIMD-aligned)
            // 3. P: Partial remainder particles at the very end (Scalar)
            const math::Range full_chunks = {range_chunks.start, range_chunks.stop - 1};
            const math::Range full_tail = {0, full_tail_end, packed_size};
            const math::Range partial_tail = {full_tail_end, tail_len};

            // Iteration Scheme for Symmetric Interaction (Upper Triangle: index of X_inta > index of X):
            // 1. [C] vs [a) Self, b) C_intra, c) C_inter, d) F, e) P]
            interact_body_vs_everything(full_chunks, full_tail, partial_tail, f);

            // 2. [F] vs [a) Self, b) F_intra, c) P]
            interact_tail_simd_vs_remainder(full_chunks, full_tail, partial_tail, full_tail_end, f);

            // 3. [P] vs [P_intra]
            interact_partial_vs_partial(full_chunks, partial_tail, tail_len, f);
        }


        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_body_vs_everything(
            const math::Range& full_chunks,
            const math::Range& full_tail,
            const math::Range& partial_tail,
            Kernel&& f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 1. all simd in full chunks vs all simd in full chunks + full tail + partial particles
            for (size_t c1 : full_chunks) {
                ctx.prefetch( c1 + 1);

                APRIL_UNROLL_LOOP_N(iter_chunks)
                for (size_t i = 0; i < chunk_size; i += packed_size) {
                    auto packed1 = ctx.packed(c1, i);
                    auto buffer1 = packed1.load_buffer();

                    // a. simd register self interaction
                    interact_symmetric_self(buffer1, packed1, f);

                    // b. intra-chunk interactions (Block i vs Blocks j where j > i inside c1)
                    for (size_t j = i + packed_size; j < chunk_size; j += packed_size) {
                        interact_block_vs_block(buffer1, ctx.packed(c1, j), f);
                    }

                    // c. inter-chunk interactions (Block i vs all Blocks in c2 where c2 > c1)
                    for (size_t c2 = c1 + 1; c2 < full_chunks.stop; ++c2) {
                        ctx.prefetch_nta( c2 + 1);
                        APRIL_UNROLL_LOOP_N(iter_chunks)
                        for (size_t j = 0; j < chunk_size; j += packed_size) {
                            interact_block_vs_block(buffer1, ctx.packed(c2, j), f);
                        }
                    }

                    // d. interaction with SIMD-aligned blocks in the tail chunk [F]
                    for (size_t j : full_tail) {
                        interact_block_vs_block(buffer1, ctx.packed(full_chunks.stop, j), f);
                    }

                    // e. interaction with partial particles in the tail chunk [P]
                    for (size_t j : partial_tail) {
                        interact_block_vs_scalar(buffer1, ctx.scalar(full_chunks.stop, j), f);
                    }

                    buffer1.update_into(packed1);
                }
            }
        }


        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_tail_simd_vs_remainder(
            const math::Range& full_chunks,
            const math::Range& full_tail,
            const math::Range& partial_tail,
            const size_t full_tail_end,
            Kernel&& f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 2. all full blocks i in tail chunk [F] vs remainder [a) Self, b) F_intra, c) P]
            for (size_t i : full_tail) {
                auto packed1 = ctx.packed(full_chunks.stop, i);
                auto buffer1 = packed1.load_buffer();

                // a. self interaction
                interact_symmetric_self(buffer1, packed1, f);

                // b. intra-chunk interactions (Block i vs Blocks j where j > i inside full tail)
                for (size_t j = i + packed_size; j < full_tail_end; j += packed_size) {
                    interact_block_vs_block(buffer1, ctx.packed(full_chunks.stop, j), f);
                }

                // c interaction with partial particles in the tail chunk [P]
                for (size_t j : partial_tail) {
                    interact_block_vs_scalar(buffer1, ctx.scalar(full_chunks.stop, j), f);
                }

                buffer1.update_into(packed1);
            }
        }


        template <typename Kernel>
        APRIL_FORCE_INLINE void interact_partial_vs_partial(
            const math::Range& full_chunks,
            const math::Range& partial_tail,
            const size_t tail_len,
            Kernel && f
        ) const {
            using namespace internal;
            BatchContext ctx(container, f);

            // 3. partial vs partial
            if (partial_tail.start != partial_tail.stop) {
                // shift lane indices to absolute indices (relative to tail chunk start)
                const auto absolute_lane_indices = packed::load_aligned(idx_arr) + static_cast<double>(partial_tail.start);

                // mask out dead sentinels
                const auto valid_tail_mask = absolute_lane_indices < static_cast<double>(tail_len);

                // load the entire tail chunk
                auto packed_chunk = ctx.packed(full_chunks.stop, partial_tail.start);
                auto chunk_data = packed_chunk.load_buffer();

                // loop over valid particles in the tail
                for (size_t i : partial_tail) {
                    auto p1 = ctx.scalar(full_chunks.stop, i);
                    auto buffer1 =  p1.broadcast();

                    // Fresh buffer to capture forces exerted on the chunk from just this p1
                    auto temp_chunk = packed_chunk.load_buffer(); // Natively zeroes WOMask

                    // scalar p1 vs entire chunk
                    auto view1 = buffer1.to_view();
                    auto view2 = temp_chunk.to_view();
                    f(view1, view2);

                    // upper triangle mask (j > i) to prevent double counting and self-interaction
                    auto current_mask = valid_tail_mask && (absolute_lane_indices > static_cast<double>(i));

                    // update p1 (scalar reduction for forces ON p1 from j > i)
                    buffer1.reduce_into(p1, current_mask);

                    // Accumulate forces ON the rest of the chunk (from p1)
                    chunk_data.accumulate(temp_chunk, current_mask);
                }
                // Safely flush the master accumulator back to memory
                chunk_data.update_into(packed_chunk);
            }
        }
    };
}
