#pragma once

#include "april/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"


namespace april::container::internal {

   	template<typename Container, typename ChunkPtr>
	struct AsymmetricChunkedBatch : SerialBatch {
		explicit AsymmetricChunkedBatch(Container & container, ChunkPtr * chunks): container(container), chunks(chunks) {}

		template<env::FieldMask Mask, typename Func>
		AP_FORCE_INLINE
   		void for_each_pair (Func && f) const {
			// skip empty range
		    if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

		    constexpr size_t stride = Container::chunk_size;

			// peel of last chunk (i.e. the tail)
		    const size_t c1_body_end = range1_chunks.stop - 1;
		    const size_t c2_body_end = range2_chunks.stop - 1;

		    // if tail is 0, it means the chunk is actually full (stride)
		    const size_t limit1_tail = (range1_tail == 0) ? stride : range1_tail;
		    const size_t limit2_tail = (range2_tail == 0) ? stride : range2_tail;

		    // body vs body with hardcoded stride x stride loop
		    for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
		        AP_PREFETCH(chunks + c1 + 1);

		        for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
		            AP_PREFETCH(chunks + c2 + 1);

		            for (size_t i = 0; i < stride; ++i) {
		                auto p1 = container.template restricted_at<Mask>(c1, i);
		                for (size_t j = 0; j < stride; ++j) {
		                    auto p2 = container.template restricted_at<Mask>(c2, j);
		                    f(p1, p2);
		                }
		            }
		        }
		    }

		    // body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
		    for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
		    	AP_PREFETCH(chunks + c1 + 1);
		        for (size_t i = 0; i < stride; ++i) {
		            auto p1 = container.template restricted_at<Mask>(c1, i);
		            for (size_t j = 0; j < limit2_tail; ++j) {
		                auto p2 = container.template restricted_at<Mask>(c2_body_end, j);
		                f(p1, p2);
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
		                f(p1, p2);
		            }
		        }
		    }

		    // tail 1 vs tail 2 (Interaction between the two last chunks)
	        for (size_t i = 0; i < limit1_tail; ++i) {
	            auto p1 = container.template restricted_at<Mask>(c1_body_end, i);
	            for (size_t j = 0; j < limit2_tail; ++j) {
	                auto p2 = container.template restricted_at<Mask>(c2_body_end, j);
	                f(p1, p2);
	            }
	        }
		}

		// Ranges represent chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
		math::Range range1_chunks; // Chunk Indices [start, end)
		size_t range1_tail{};  // Number of valid items in the last chunk of range1 (0 = Full)

		math::Range range2_chunks;
		size_t range2_tail{};
	private:
		Container & container;
		ChunkPtr * AP_RESTRICT const chunks = container.ptr_chunks;
	};


	template<typename Container, typename ChunkPtr>
	struct SymmetricChunkedBatch : SerialBatch {
		explicit SymmetricChunkedBatch(Container & container, ChunkPtr * chunks) : container(container), chunks(chunks) {}

		template<env::FieldMask Mask, typename Func>
	    AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
	        if (range_chunks.start == range_chunks.stop) return;
	        constexpr size_t stride = Container::chunk_size;

			// peel of last chunk (i.e. the tail)
	        const size_t c_body_end = range_chunks.stop - 1;
	        const size_t limit_tail = (range_tail == 0) ? stride : range_tail;

	        // body (iterate c1 up to the last full chunk)
	        for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
	            AP_PREFETCH(chunks + c1 + 1);

	            // chunk self interaction
	            for (size_t i = 0; i < stride; ++i) {
	                auto p1 = container.template restricted_at<Mask>(c1, i);
	                for (size_t j = i + 1; j < stride; ++j) {
	                     auto p2 = container.template restricted_at<Mask>(c1, j);
	                     f(p1, p2);
	                }
	            }

	        	// interaction with all other loops
	            for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
	                AP_PREFETCH(chunks + c2 + 1);
	                for (size_t i = 0; i < stride; ++i) {
	                    auto p1 = container.template restricted_at<Mask>(c1, i);
	                    for (size_t j = 0; j < stride; ++j) {
	                        auto p2 = container.template restricted_at<Mask>(c2, j);
	                        f(p1, p2);
	                    }
	                }
	            }
	        }

	        // body vs tail (every body chunk with tail chunk)
	        for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
	            AP_PREFETCH(chunks + c1 + 1);
	            for (size_t i = 0; i < stride; ++i) {
	                auto p1 = container.template restricted_at<Mask>(c1, i);
	                for (size_t j = 0; j < limit_tail; ++j) {
	                     auto p2 = container.template restricted_at<Mask>(c_body_end, j);
	                     f(p1, p2);
	                }
	            }
	        }

	        // tail (interact tail chunk with itself)
            for (size_t i = 0; i < limit_tail; ++i) {
                auto p1 = container.template restricted_at<Mask>(c_body_end, i);
                for (size_t j = i + 1; j < limit_tail; ++j) {
                     auto p2 = container.template restricted_at<Mask>(c_body_end, j);
                     f(p1, p2);
                }
            }
	    }

		// Range represents chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
		math::Range  range_chunks;  // Chunk Indices [start, end)
		size_t range_tail{};  // Number of valid items in the last chunk of range1 (0 = Full)
	private:
		Container & container;
		ChunkPtr * AP_RESTRICT const chunks = container.ptr_chunks;
	};

}

