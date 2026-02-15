#pragma once
#include "april/base/types.hpp"
#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"

#include "april/simd/packed.hpp"
#include "april/math/range.hpp"

namespace april::container {
    template<typename Container, typename ChunkPtr>
    struct AsymmetricChunkedSimdBatch : BatchBase<ParallelPolicy::None, UpdatePolicy::Serial, ComputePolicy::Vector> {

        // We need to know the width of our SIMD registers
        static constexpr size_t stride = packed::size();
        static constexpr size_t ChunkSize = Container::chunk_size;
        static constexpr size_t Iterations = ChunkSize / stride;
        static_assert(ChunkSize % stride == 0, "Chunk Size must be a multiple of SIMD Width");
        static_assert(ChunkSize == stride, "For now: Chunk size == SIMD WIDTH");

        explicit AsymmetricChunkedSimdBatch(Container & container, ChunkPtr * chunks)
            : container(container), chunks(chunks) {}

        template<env::Field Mask, typename Func>
        AP_FORCE_INLINE
        void for_each_pair(Func && kernel) const {
            // Skip empty range
            if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

            // Peel off last chunk (Tail)
            const size_t c1_body_end = range1_chunks.stop - 1;
            const size_t c2_body_end = range2_chunks.stop - 1;

            // Calculate actual tail counts (0 means full chunk)
            const size_t limit1_tail = (range1_tail == 0) ? ChunkSize : range1_tail;
            const size_t limit2_tail = (range2_tail == 0) ? ChunkSize : range2_tail;

            for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
                AP_PREFETCH(chunks + c1 + 1);
                auto packed1 = container.template restricted_at_packed<Mask>(c1, 0);
                auto buffer1 = packed1.load_buffer();

                for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
                    AP_PREFETCH(chunks + c2 + 1);
                    auto packed2 = container.template restricted_at_packed<Mask>(c2, 0);
                    auto buffer2 = packed2.load_buffer();

                    for (size_t k = 0; k < stride; k++) {
                    	kernel(buffer1, buffer2);
                        buffer2.rotate_right();
                    }
                    packed2.force += buffer2.force;
                }

                packed1.force += buffer1.force;
            }


		    // body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
		    for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
		    	AP_PREFETCH(chunks + c1 + 1);
		        for (size_t i = 0; i < stride; ++i) {
		            auto p1 = container.template restricted_at<Mask>(c1, i);
		            for (size_t j = 0; j < limit2_tail; ++j) {
		                auto p2 = container.template restricted_at<Mask>(c2_body_end, j);
		            	kernel(p1, p2);
		            }
		        }
		    }

			// body 2 vs tail 1 (iterate full chunks of R2 against the single partial chunk of R1
		    for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
		    	AP_PREFETCH(chunks + c2 + 1);
		        for (size_t i = 0; i < limit1_tail; ++i) {
		            auto p1 = container.template restricted_at<Mask>(c1_body_end, i);
		            for (size_t j = 0; j < stride; ++j) {
		                auto p2 = container.template restricted_at<Mask>(c2, j);
		            	kernel(p1, p2);
		            }
		        }
		    }

		    // tail 1 vs tail 2 (Interaction between the two last chunks)
	        for (size_t i = 0; i < limit1_tail; ++i) {
	            auto p1 = container.template restricted_at<Mask>(c1_body_end, i);
	            for (size_t j = 0; j < limit2_tail; ++j) {
	                auto p2 = container.template restricted_at<Mask>(c2_body_end, j);
	            	kernel(p1, p2);
	            }
	        }
        }

        // Members matching the scalar batch
        math::Range range1_chunks;
        size_t range1_tail{};

        math::Range range2_chunks;
        size_t range2_tail{};

    private:
        Container & container;
        ChunkPtr * AP_RESTRICT const chunks;
    };

}


