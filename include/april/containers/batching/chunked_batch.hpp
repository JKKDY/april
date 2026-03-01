#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"


namespace april::container::batching {
    //-----------------
    // ASYMMETRIC BATCH
    //-----------------
    template <typename Container, typename ChunkPtr>
    struct AsymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
        explicit
        AsymmetricChunkedBatch(Container& container, ChunkPtr* chunks) : container(container), chunks(chunks) {
            for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
        }

        template <ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            // skip empty range
            if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

            if constexpr (static_cast<bool>(E & exec::internal::ExecutionMode::Vector)) {
                for_each_pair_packed<P>(f);
            }
            else {
                for_each_pair_scalar<P>(f);
            }
        }

        // range1_chunks/range2_chunks represent chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
        math::Range range1_chunks; // Chunk Indices [start, end)
        math::Range range2_chunks;
        size_t range1_tail{}; // Number of valid items in the last chunk of range1 (0 = Full)
        size_t range2_tail{};

    private:
        Container& container;
        ChunkPtr* AP_RESTRICT const chunks = container.ptr_chunks;

        static constexpr size_t chunk_size = Container::chunk_size;
        static constexpr size_t packed_size = packed::size();
        static constexpr size_t iter_chunks = chunk_size / packed_size;

        alignas(64) double idx_arr[packed_size];


        //----------------
        // SCALAR ITERATOR
        //----------------
        template <ParallelPolicy P, exec::IsKernel Kernel>
        void for_each_pair_scalar(Kernel && f) const {
            using K = std::remove_cvref_t<Kernel>;
            // peel of last chunk (i.e. the tail)
            const size_t c1_body_end = range1_chunks.stop - 1;
            const size_t c2_body_end = range2_chunks.stop - 1;

            // if tail is 0, it means the chunk is actually full (stride)
            const size_t limit1_tail = (range1_tail == 0) ? chunk_size : range1_tail;
            const size_t limit2_tail = (range2_tail == 0) ? chunk_size : range2_tail;

            // body vs body (
            for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
                AP_PREFETCH(chunks + c1 + 1);

                for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
                    AP_PREFETCH(chunks + c2 + 1);

                    for (size_t i = 0; i < chunk_size; ++i) {
                        auto p1 = container.template at<K::access>(c1, i);
                        for (size_t j = 0; j < chunk_size; ++j) {
                            auto p2 = container.template at<K::access>(c2, j);
                            f(p1, p2);
                        }
                    }
                }
            }

            // body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
            for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
                AP_PREFETCH(chunks + c1 + 1);
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = container.template at<K::access>(c1, i);
                    for (size_t j = 0; j < limit2_tail; ++j) {
                        auto p2 = container.template at<K::access>(c2_body_end, j);
                        f(p1, p2);
                    }
                }
            }

            // body 2 vs tail 1 (iterate full chunks of R2 against the single partial chunk of R1
            for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
                AP_PREFETCH(chunks + c2 + 1);
                for (size_t i = 0; i < limit1_tail; ++i) {
                    auto p1 = container.template at<K::access>(c1_body_end, i);
                    for (size_t j = 0; j < chunk_size; ++j) {
                        auto p2 = container.template at<K::access>(c2, j);
                        f(p1, p2);
                    }
                }
            }

            // tail 1 vs tail 2 (Interaction between the two last chunks)
            for (size_t i = 0; i < limit1_tail; ++i) {
                auto p1 = container.template at<K::access>(c1_body_end, i);
                for (size_t j = 0; j < limit2_tail; ++j) {
                    auto p2 = container.template at<K::access>(c2_body_end, j);
                    f(p1, p2);
                }
            }
        }


        //----------------
        // PACKED ITERATOR
        //----------------
        template <ParallelPolicy P, exec::IsKernel Kernel>
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
        AP_FORCE_INLINE void interact_body_vs_body(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& full_tail1, const math::Range& full_tail2,
            Kernel && f) const {

            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 1. full SIMD blocks in range 1 (Chunks) vs full SIMD blocks in range 2 (Chunks + Tail)
            for (size_t c1 : full_chunks1) {
                AP_PREFETCH(chunks + c1 + 1);
                AP_UNROLL_LOOP_N(iter_chunks)
                for (size_t i = 0; i < chunk_size; i += packed_size) {
                    auto packed1 = container.template at_packed<access_fields>(c1, i);
                    auto buffer1 = packed1.load_buffer(); // load simd block in range 1

                    // a. sweep buffer1 across all full SIMD blocks in Range 2 body chunks [C2]
                    for (size_t c2 : full_chunks2) {
                        AP_PREFETCH(chunks + c2 + 1);
                        AP_UNROLL_LOOP_N(iter_chunks)
                        for (size_t j = 0; j < chunk_size; j += packed_size) {
                            auto packed2 = container.template at_packed<access_fields>(c2, j);
                            auto buffer2 = packed2.load_buffer();

                            AP_UNROLL_LOOP_N(packed_size)
                            for (size_t k = 0; k < packed_size; k++) {
                                f(buffer1, buffer2);
                                buffer2.rotate_right();
                            }

                            packed2.force = buffer2.force;
                        }
                    }

                    // b. sweep buffer1 across SIMD-aligned blocks within the Range 2 tail chunk [F2]
                    for (size_t t2 : full_tail2) {
                        auto packed2 = container.template at_packed<access_fields>(full_chunks2.stop, t2);
                        auto buffer2 = packed2.load_buffer();

                        AP_UNROLL_LOOP_N(packed_size)
                        for (size_t k = 0; k < packed_size; k++) {
                            f(buffer1, buffer2);
                            buffer2.rotate_right();
                        }

                        packed2.force = buffer2.force;
                    }

                    packed1.force = buffer1.force; // final register flush for buffer1
                }
            }

            // 2. SIMD-aligned blocks in the Range 1 tail chunk vs Full Range 2 (Body + Tail)
            for (size_t t1 : full_tail1) {
                auto packed1 = container.template at_packed<access_fields>(full_chunks1.stop, t1);
                auto buffer1 = packed1.load_buffer();

                // a. Interaction with all full chunks in the Range 2 body [C2]
                for (size_t c2 : full_chunks2) {
                    AP_PREFETCH(chunks + c2 + 1);
                    AP_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size) {
                        auto packed2 = container.template at_packed<access_fields>(c2, j);
                        auto buffer2 = packed2.load_buffer();

                        AP_UNROLL_LOOP_N(packed_size)
                        for (size_t k = 0; k < packed_size; k++) {
                            f(buffer1, buffer2);
                            buffer2.rotate_right();
                        }

                        packed2.force = buffer2.force;
                    }
                }

                // b. Interaction with SIMD-aligned blocks within the Range 2 tail chunk [F2]
                for (size_t t2 : full_tail2) {
                    auto packed2 = container.template at_packed<access_fields>(full_chunks2.stop, t2);
                    auto buffer2 = packed2.load_buffer();

                    AP_UNROLL_LOOP_N(packed_size)
                    for (size_t k = 0; k < packed_size; k++) {
                        f(buffer1, buffer2);
                        buffer2.rotate_right();
                    }

                    packed2.force = buffer2.force;
                }
                packed1.force = buffer1.force;
            }
        }


        template <typename Kernel>
        AP_FORCE_INLINE void interact_tails_vs_body(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& full_tail1, const math::Range& full_tail2,
            const math::Range& partial_tail1, const math::Range& partial_tail2,
            Kernel && f) const {

            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 3. Range 2 Partial Tail [P2] vs Range 1 SIMD blocks [C1 + F1]
            for (size_t i : partial_tail2) {
                auto p2 = container.template at<access_fields>(full_chunks2.stop, i);
                auto buffer2 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p2);
                buffer2.force = {0, 0, 0};

                // a. Interaction with full chunks in Range 1 body
                for (size_t c1 : full_chunks1) {
                    AP_PREFETCH(chunks + c1 + 1);
                    AP_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size) {
                        auto packed1 = container.template at_packed<access_fields>(c1, j);
                        auto buffer1 = packed1.load_buffer();
                        buffer1.force = {0, 0, 0};

                        f(buffer1, buffer2);

                        packed1.force += buffer1.force;
                    }
                }

                // b. Interaction with SIMD-aligned blocks in Range 1 tail chunk
                for (size_t t1 : full_tail1) {
                    auto packed1 = container.template at_packed<access_fields>(full_chunks1.stop, t1);
                    auto buffer1 = packed1.load_buffer();
                    buffer1.force = {0, 0, 0};

                    f(buffer1, buffer2);

                    packed1.force += buffer1.force;
                }

                p2.force.x += buffer2.force.x.reduce_add();
                p2.force.y += buffer2.force.y.reduce_add();
                p2.force.z += buffer2.force.z.reduce_add();
            }

            // Partial Tail 1 vs Full Range 2 (chunks + full tail)
            for (size_t i : partial_tail1) {
                auto p1 = container.template at<access_fields>(full_chunks1.stop, i);
                auto buffer1 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p1);
                buffer1.force = {0, 0, 0};

                // a. Interaction with full chunks in Range 2 body
                for (size_t c2 : full_chunks2) {
                    AP_PREFETCH(chunks + c2 + 1);
                    AP_UNROLL_LOOP_N(iter_chunks)
                    for (size_t j = 0; j < chunk_size; j += packed_size) {
                        auto packed2 = container.template at_packed<access_fields>(c2, j);
                        auto buffer2 = packed2.load_buffer();
                        buffer2.force = {0, 0, 0};

                        f(buffer1, buffer2);

                        packed2.force += buffer2.force;
                    }
                }

                // b. Interaction with SIMD-aligned blocks in Range 2 tail chunk
                for (size_t t2 : full_tail2) {
                    auto packed2 = container.template at_packed<access_fields>(full_chunks2.stop, t2);
                    auto buffer2 = packed2.load_buffer();
                    buffer2.force = {0, 0, 0};

                    f(buffer1, buffer2);

                    packed2.force += buffer2.force;
                }

                p1.force.x += buffer1.force.x.reduce_add();
                p1.force.y += buffer1.force.y.reduce_add();
                p1.force.z += buffer1.force.z.reduce_add();
            }
        }

        template <typename Kernel>
        AP_FORCE_INLINE void interact_tails_vs_tails(
            const math::Range& full_chunks1, const math::Range& full_chunks2,
            const math::Range& partial_tail1, const math::Range& partial_tail2,
            Kernel && f) const {
            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 5. partial tail vs partial tail
            if (partial_tail1.start != partial_tail1.stop && partial_tail2.start != partial_tail2.stop) {
                // initialize mask for the valid lanes in partial_tail1
                const auto lane_indices = packed::load_aligned(idx_arr);

                const double valid_lanes = static_cast<double>(partial_tail1.stop - partial_tail1.start);
                const auto mask = lane_indices < valid_lanes;
                const packed null = 0.0;

                // load the single SIMD block containing partial_tail1
                auto packed1 = container.template at_packed<access_fields>(full_chunks1.stop, partial_tail1.start);
                auto buffer1 = packed1.load_buffer();
                buffer1.force = {0, 0, 0};

                // loop over the individual particles in partial_tail2
                for (size_t i = partial_tail2.start; i < partial_tail2.stop; ++i) {
                    auto p2 = container.template at<access_fields>(full_chunks2.stop, i);
                    auto buffer2 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p2);
                    buffer2.force = {0, 0, 0};

                    f(buffer1, buffer2);

                    p2.force += vec3{
                        select(mask, buffer2.force.x, null).reduce_add(),
                        select(mask, buffer2.force.y, null).reduce_add(),
                        select(mask, buffer2.force.z, null).reduce_add()
                    };
                }

                // write back valid lanes to memory for packed1
                packed1.force += pvec3{
                    select(mask, buffer1.force.x, null),
                    select(mask, buffer1.force.y, null),
                    select(mask, buffer1.force.z, null)
                };
            }
        }
    };




    //================
    //----------------
    // SYMMETRIC BATCH
    //----------------
    //================
    template <typename Container, typename ChunkPtr>
    struct SymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
        explicit SymmetricChunkedBatch(Container& container, ChunkPtr* chunks) : container(container), chunks(chunks) {
            for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
        }

        template <ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair(Func&& f) const {
            if (range_chunks.start == range_chunks.stop) return;

            if constexpr (static_cast<bool>(E & exec::internal::ExecutionMode::Vector)) {
                for_each_pair_packed<P>(f);
            }
            else {
                for_each_pair_scalar<P>(f);
            }
        }

        // Range represents chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
        math::Range range_chunks; // Chunk Indices [start, end)
        size_t range_tail{}; // Number of valid items in the last chunk of range1 (0 = Full)
    private:
        Container& container;
        ChunkPtr* AP_RESTRICT const chunks = container.ptr_chunks;

        static constexpr size_t chunk_size = Container::chunk_size;
        static constexpr size_t packed_size = packed::size();
        static constexpr size_t iter_chunks = chunk_size / packed_size;

        alignas(64) double idx_arr[packed_size];


        //----------------
        // SCALAR ITERATOR
        //----------------
        template <ParallelPolicy P, typename Kernel>
        void for_each_pair_scalar(Kernel && f) const {
            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            const size_t c_body_end = range_chunks.stop - 1;
            const size_t limit_tail = (range_tail == 0) ? packed_size : range_tail;

            for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
                AP_PREFETCH(chunks + c1 + 1);

                // chunk self interaction
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = container.template at<access_fields>(c1, i);
                    for (size_t j = i + 1; j < chunk_size; ++j) {
                        auto p2 = container.template at<access_fields>(c1, j);
                        f(p1, p2);
                    }
                }

                // interaction with all other chunks
                for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
                    AP_PREFETCH(chunks + c2 + 1);
                    for (size_t i = 0; i < chunk_size; ++i) {
                        auto p1 = container.template at<access_fields>(c1, i);
                        for (size_t j = 0; j < chunk_size; ++j) {
                            auto p2 = container.template at<access_fields>(c2, j);
                            f(p1, p2);
                        }
                    }
                }
            }

            // body vs tail (every body chunk with tail chunk)
            for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
                AP_PREFETCH(chunks + c1 + 1);
                for (size_t i = 0; i < chunk_size; ++i) {
                    auto p1 = container.template at<access_fields>(c1, i);
                    for (size_t j = 0; j < limit_tail; ++j) {
                        auto p2 = container.template at<access_fields>(c_body_end, j);
                        f(p1, p2);
                    }
                }
            }

            // tail (interact tail chunk with itself)
            for (size_t i = 0; i < limit_tail; ++i) {
                auto p1 = container.template at<access_fields>(c_body_end, i);
                for (size_t j = i + 1; j < limit_tail; ++j) {
                    auto p2 = container.template at<access_fields>(c_body_end, j);
                    f(p1, p2);
                }
            }
        }


        //----------------
        // PACKED ITERATOR
        //----------------
        template <ParallelPolicy P>
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
        AP_FORCE_INLINE void interact_body_vs_everything(
            const math::Range& full_chunks,
            const math::Range& full_tail,
            const math::Range& partial_tail,
            Kernel&& f) const {
            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 1. all simd in full chunks vs all simd in full chunks + full tail + partial particles
            for (size_t c1 : full_chunks) {
                AP_PREFETCH(chunks + c1 + 1);

                AP_UNROLL_LOOP_N(iter_chunks)
                for (size_t i = 0; i < chunk_size; i += packed_size) {
                    auto packed1 = container.template at_packed<access_fields>(c1, i);
                    auto buffer1 = packed1.load_buffer();

                    // a. simd register self interaction
                    {
                        auto buffer2 = packed1.load_buffer();
                        buffer2.force = {0, 0, 0};

                        AP_UNROLL_LOOP()
                        for (size_t k = 0; k < packed_size / 2 - 1; k++) {
                            buffer2.rotate_right();
                            f(buffer1, buffer2);
                        }

                        buffer1.force.x += buffer2.force.x.template rotate_left<packed_size / 2 - 1>();
                        buffer1.force.y += buffer2.force.y.template rotate_left<packed_size / 2 - 1>();
                        buffer1.force.z += buffer2.force.z.template rotate_left<packed_size / 2 - 1>();

                        packed1.force = buffer1.force;

                        buffer2.rotate_right();
                        f(buffer1, buffer2);
                    }

                    // b. intra-chunk interactions (Block i vs Blocks j where j > i inside c1)
                    for (size_t j = i + packed_size; j < chunk_size; j += packed_size) {
                        auto packed2 = container.template at_packed<access_fields>(c1, j);
                        auto buffer2 = packed2.load_buffer();

                        AP_UNROLL_LOOP_N(packed_size)
                        for (size_t k = 0; k < packed_size; k++) {
                            f(buffer1, buffer2);
                            buffer2.template rotate_right<1>();
                        }
                        packed2.force = buffer2.force;
                    }

                    // c. inter-chunk interactions (Block i vs all Blocks in c2 where c2 > c1)
                    for (size_t c2 = c1 + 1; c2 < full_chunks.stop; ++c2) {
                        AP_PREFETCH(chunks + c2 + 1);
                        AP_UNROLL_LOOP_N(iter_chunks)
                        for (size_t j = 0; j < chunk_size; j += packed_size) {
                            auto packed2 = container.template at_packed<access_fields>(c2, j);
                            auto buffer2 = packed2.load_buffer();

                            AP_UNROLL_LOOP_N(packed_size)
                            for (size_t k = 0; k < packed_size; k++) {
                                f(buffer1, buffer2);
                                buffer2.template rotate_right<1>();
                            }
                            packed2.force = buffer2.force;
                        }
                    }

                    // d. interaction with SIMD-aligned blocks in the tail chunk [F]
                    for (size_t j : full_tail) {
                        auto packed2 = container.template at_packed<access_fields>(full_chunks.stop, j);
                        auto buffer2 = packed2.load_buffer();

                        AP_UNROLL_LOOP_N(packed_size)
                        for (size_t k = 0; k < packed_size; k++) {
                            f(buffer1, buffer2);
                            buffer2.template rotate_right<1>();
                        }

                        packed2.force = buffer2.force;
                    }

                    // e. interaction with partial particles in the tail chunk [P]
                    for (size_t j : partial_tail) {
                        auto p2 = container.template at<access_fields>(full_chunks.stop, j);
                        auto buffer2 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p2);
                        buffer2.force = {0, 0, 0};

                        f(buffer1, buffer2);

                        p2.force.x += buffer2.force.x.reduce_add();
                        p2.force.y += buffer2.force.y.reduce_add();
                        p2.force.z += buffer2.force.z.reduce_add();
                    }

                    packed1.force = buffer1.force;
                }
            }
        }


        template <typename Kernel>
        AP_FORCE_INLINE void interact_tail_simd_vs_remainder(
            const math::Range& full_chunks,
            const math::Range& full_tail,
            const math::Range& partial_tail,
            const size_t full_tail_end,
            Kernel&& f) const {

            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 2. all full blocks i in tail chunk [F] vs remainder [a) Self, b) F_intra, c) P]
            for (size_t i : full_tail) {
                // a. self interaction
                auto packed1 = container.template at_packed<access_fields>(full_chunks.stop, i);
                auto buffer1 = packed1.load_buffer();
                {
                    auto buffer2 = packed1.load_buffer();
                    buffer2.force = {0, 0, 0};

                    AP_UNROLL_LOOP()
                    for (size_t k = 0; k < packed_size / 2 - 1; k++) {
                        buffer2.rotate_right();
                        f(buffer1, buffer2);
                    }

                    buffer1.force.x += buffer2.force.x.template rotate_left<packed_size / 2 - 1>();
                    buffer1.force.y += buffer2.force.y.template rotate_left<packed_size / 2 - 1>();
                    buffer1.force.z += buffer2.force.z.template rotate_left<packed_size / 2 - 1>();

                    buffer2.rotate_right();
                    f(buffer1, buffer2);
                }

                // b. intra-chunk interactions (Block i vs Blocks j where j > i inside full tail)
                for (size_t j = i + packed_size; j < full_tail_end; j += packed_size) {
                    auto packed2 = container.template at_packed<access_fields>(full_chunks.stop, j);
                    auto buffer2 = packed2.load_buffer();

                    AP_UNROLL_LOOP_N(packed_size)
                    for (size_t k = 0; k < packed_size; k++) {
                        f(buffer1, buffer2);
                        buffer2.rotate_right();
                    }
                    packed2.force = buffer2.force;
                }

                // c interaction with partial particles in the tail chunk [P]
                for (size_t j : partial_tail) {
                    auto p2 = container.template at<access_fields>(full_chunks.stop, j);
                    auto buffer2 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p2);
                    buffer2.force = {0, 0, 0};

                    f(buffer1, buffer2);

                    p2.force.x += buffer2.force.x.reduce_add();
                    p2.force.y += buffer2.force.y.reduce_add();
                    p2.force.z += buffer2.force.z.reduce_add();
                }

                packed1.force = buffer1.force;
            }
        }

        template <typename Kernel>
        AP_FORCE_INLINE void interact_partial_vs_partial(
            const math::Range& full_chunks,
            const math::Range& partial_tail,
            const size_t tail_len,
            Kernel && f) const {

            using K = std::remove_cvref_t<Kernel>;
            constexpr ParticleField access_fields = K::access;

            // 3. partial vs partial
            if (partial_tail.start != partial_tail.stop) {
                // shift lane indices to absolute indices (relative to tail chunk start)
                const auto absolute_lane_indices = packed::load_aligned(idx_arr) + static_cast<double>(partial_tail.
                    start);

                // mask out dead sentinels
                const auto valid_tail_mask = absolute_lane_indices < static_cast<double>(tail_len);
                const packed null = 0.0;

                // load the entire tail chunk
                auto packed_chunk = container.template at_packed<access_fields>(full_chunks.stop, partial_tail.start);
                auto chunk_data = packed_chunk.load_buffer();

                // Accumulator for forces exerted by the p1 scalars onto the rest of the chunk
                pvec3 accum_force = {0, 0, 0};

                // loop over valid particles in the tail
                for (size_t i : partial_tail) {
                    auto p1 = container.template at<access_fields>(full_chunks.stop, i);
                    auto buffer1 = particle::internal::PackedParticleBuffer<access_fields>::broadcast(p1);
                    buffer1.force = {0, 0, 0};

                    // Fresh buffer to capture forces exerted on the chunk from just this p1
                    auto temp_chunk = chunk_data;
                    temp_chunk.force = {0, 0, 0};

                    // scalar p1 vs entire chunk
                    f(buffer1, temp_chunk);

                    // FIX: Compare against absolute_lane_indices
                    auto current_mask = valid_tail_mask && (absolute_lane_indices > static_cast<double>(i));

                    // update p1 (scalar reduction for forces ON p1 from j > i)
                    p1.force += vec3{
                        select(current_mask, buffer1.force.x, null).reduce_add(),
                        select(current_mask, buffer1.force.y, null).reduce_add(),
                        select(current_mask, buffer1.force.z, null).reduce_add()
                    };

                    // Accumulate forces ON the rest of the chunk (from p1)
                    accum_force.x += select(current_mask, temp_chunk.force.x, null);
                    accum_force.y += select(current_mask, temp_chunk.force.y, null);
                    accum_force.z += select(current_mask, temp_chunk.force.z, null);
                }

                packed_chunk.force += accum_force;
            }
        }
    };
}
