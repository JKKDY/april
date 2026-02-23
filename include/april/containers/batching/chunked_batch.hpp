#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"


namespace april::container::internal {

   	template<typename Container, typename ChunkPtr>
	struct AsymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
		explicit AsymmetricChunkedBatch(Container & container, ChunkPtr * chunks): container(container), chunks(chunks) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
		AP_FORCE_INLINE
   		void for_each_pair (Func && f) const {
			// skip empty range
		    if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

		    constexpr size_t simd_width = Container::chunk_size;
			constexpr bool vectorize = static_cast<bool>(E & exec::internal::ExecutionMode::Vector);

			static_assert(packed::size() == Container::chunk_size);

			// peel of last chunk (i.e. the tail)
		    const size_t c1_body_end = range1_chunks.stop - 1;
		    const size_t c2_body_end = range2_chunks.stop - 1;

		    // if tail is 0, it means the chunk is actually full (stride)
		    const size_t limit1_tail = (range1_tail == 0) ? simd_width : range1_tail;
		    const size_t limit2_tail = (range2_tail == 0) ? simd_width : range2_tail;

		    // body vs body with hardcoded stride x stride loop
			if constexpr (vectorize) {
				for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
					AP_PREFETCH(chunks + c1 + 1);
					auto packed1 = container.template at_packed<Mask>(c1, 0);
					auto buffer1 = packed1.load_buffer();

					for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
						AP_PREFETCH(chunks + c2 + 1);
						auto packed2 = container.template at_packed<Mask>(c2, 0);
						auto buffer2 = packed2.load_buffer();

						AP_UNROLL_LOOP_N(simd_width)
						for (size_t k = 0; k < simd_width; k++) {
							f(buffer1, buffer2);
							buffer2.rotate_right();
						}

						packed2.force = buffer2.force;
					}

					packed1.force = buffer1.force;
				}
			} else {
			    for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
			        AP_PREFETCH(chunks + c1 + 1);

			        for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
			            AP_PREFETCH(chunks + c2 + 1);

			            for (size_t i = 0; i < simd_width; ++i) {
			                auto p1 = container.template at<Mask>(c1, i);
			                for (size_t j = 0; j < simd_width; ++j) {
			                    auto p2 = container.template at<Mask>(c2, j);
			                    f(p1, p2);
			                }
			            }
			        }
			    }
			}

			// body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
			if constexpr (vectorize) {
				for (size_t i = 0; i < limit2_tail; ++i) {
					auto p2 = container.template at<Mask>(c2_body_end, i);
					auto buffer2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
					buffer2.force = {0,0,0};

					for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
						AP_PREFETCH(chunks + c1 + 1);

						auto packed1 = container.template at_packed<Mask>(c1, 0);
						auto buffer1 = packed1.load_buffer();
						buffer1.force = {0,0,0};

						f(buffer1, buffer2);

						packed1.force += buffer1.force;
					}

					p2.force.x += buffer2.force.x.reduce_add();
					p2.force.y += buffer2.force.y.reduce_add();
					p2.force.z += buffer2.force.z.reduce_add();
				}
			} else {
				for (size_t c1 = range1_chunks.start; c1 < c1_body_end; ++c1) {
					AP_PREFETCH(chunks + c1 + 1);
					for (size_t i = 0; i < simd_width; ++i) {
						auto p1 = container.template at<Mask>(c1, i);
						for (size_t j = 0; j < limit2_tail; ++j) {
							auto p2 = container.template at<Mask>(c2_body_end, j);
							f(p1, p2);
						}
					}
				}
			}

			// body 2 vs tail 1 (iterate full chunks of R2 against the single partial chunk of R1
			if constexpr (vectorize) {
				for (size_t i = 0; i < limit1_tail; ++i) {
					auto p1 = container.template at<Mask>(c1_body_end, i);
					auto buffer1 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p1);
					buffer1.force = {0,0,0};

					for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
						AP_PREFETCH(chunks + c2 + 1);

						auto packed2 = container.template at_packed<Mask>(c2, 0);
						auto buffer2 = packed2.load_buffer();
						buffer2.force = {0,0,0};

						f(buffer1, buffer2);

						packed2.force += buffer2.force;
					}

					p1.force.x += buffer1.force.x.reduce_add();
					p1.force.y += buffer1.force.y.reduce_add();
					p1.force.z += buffer1.force.z.reduce_add();
				}
			} else {
				for (size_t c2 = range2_chunks.start; c2 < c2_body_end; ++c2) {
					AP_PREFETCH(chunks + c2 + 1);
					for (size_t i = 0; i < limit1_tail; ++i) {
						auto p1 = container.template at<Mask>(c1_body_end, i);
						for (size_t j = 0; j < simd_width; ++j) {
							auto p2 = container.template at<Mask>(c2, j);
							f(p1, p2);
						}
					}
				}
			}

			if constexpr (static_cast<bool>(E & exec::internal::ExecutionMode::Vector)) {
				// initialize mask
				alignas(64) double idx_arr[simd_width];
				for (size_t k = 0; k < simd_width; ++k) idx_arr[k] = static_cast<double>(k);
				auto lane_indices = packed::load_aligned(idx_arr);
				auto mask = lane_indices < static_cast<double>(limit1_tail);

				const packed null = 0.0;

				// load first chunk
				auto packed1 = container.template at_packed<Mask>(c1_body_end, 0);
				auto buffer1 = packed1.load_buffer();
				buffer1.force = {0, 0, 0};

				// loop over second chunk till cutoff
				for (size_t i = 0; i < limit2_tail; ++i) {
					auto p2 = container.template at<Mask>(c2_body_end, i);
					auto buffer2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
					buffer2.force = {0, 0, 0};

					// scalar p2 vs entire Chunk 1
					f(buffer1, buffer2);

					// update p2
					p2.force += vec3 {
						select(mask, buffer2.force.x, null).reduce_add(),
						select(mask, buffer2.force.y, null).reduce_add(),
						select(mask, buffer2.force.z, null).reduce_add()
					 };
				}

				// update chunk 1
				packed1.force +=pvec3 {
					select(mask, buffer1.force.x, null),
					select(mask, buffer1.force.y, null),
					select(mask, buffer1.force.z, null)
				};
			} else {
				// tail 1 vs tail 2 (Interaction between the two last chunks)
				for (size_t i = 0; i < limit1_tail; ++i) {
					auto p1 = container.template at<Mask>(c1_body_end, i);
					for (size_t j = 0; j < limit2_tail; ++j) {
						auto p2 = container.template at<Mask>(c2_body_end, j);
						f(p1, p2);
					}
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
	struct SymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
		explicit SymmetricChunkedBatch(Container & container, ChunkPtr * chunks) : container(container), chunks(chunks) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
	    AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
	        if (range_chunks.start == range_chunks.stop) return;

	        constexpr size_t simd_width = Container::chunk_size;
			constexpr bool vectorize = static_cast<bool>(E & exec::internal::ExecutionMode::Vector);

			// peel of last chunk (i.e. the tail)
	        const size_t c_body_end = range_chunks.stop - 1;
	        const size_t limit_tail = (range_tail == 0) ? simd_width : range_tail;

	        // body (iterate c1 up to the last full chunk)
			if constexpr (vectorize) {
				for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
					AP_PREFETCH(chunks + c1 + 1);

					// chunk self interaction
					{
						auto packed1 = container.template at_packed<Mask>(c1, 0);
						auto buffer1 = packed1.load_buffer();
						buffer1.force = {0,0,0};
						auto buffer2 = packed1.load_buffer();
						buffer2.force = {0,0,0};

						AP_UNROLL_LOOP()
						for (size_t i = 0; i < simd_width/2 - 1; i++) {
							buffer2.rotate_right();
							f(buffer1, buffer2);
						}

						packed1.force.x += buffer2.force.x.template rotate_left<simd_width/2 - 1>();
						packed1.force.y += buffer2.force.y.template rotate_left<simd_width/2 - 1>();
						packed1.force.z += buffer2.force.z.template rotate_left<simd_width/2 - 1>();

						buffer2.rotate_right();
						f(buffer1, buffer2);
					}

					auto packed1 = container.template at_packed<Mask>(c1, 0);
					auto buffer1 = packed1.load_buffer();

					for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
						AP_PREFETCH(chunks + c2 + 1);
						auto packed2 = container.template at_packed<Mask>(c2, 0);
						auto buffer2 = packed2.load_buffer();

						AP_UNROLL_LOOP_N(simd_width)
						for (size_t k = 0; k < simd_width; k++) {
							f(buffer1, buffer2);
							buffer2.template rotate_right<1>();
						}
						packed2.force = buffer2.force;
					}
					packed1.force = buffer1.force;
				}

			} else {
		        for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
		            AP_PREFETCH(chunks + c1 + 1);

		            // chunk self interaction
		            for (size_t i = 0; i < simd_width; ++i) {
		                auto p1 = container.template at<Mask>(c1, i);
		                for (size_t j = i + 1; j < simd_width; ++j) {
		                     auto p2 = container.template at<Mask>(c1, j);
		                     f(p1, p2);
		                }
		            }

	        		// interaction with all other chunks
		            for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
		                AP_PREFETCH(chunks + c2 + 1);
		                for (size_t i = 0; i < simd_width; ++i) {
		                    auto p1 = container.template at<Mask>(c1, i);
		                    for (size_t j = 0; j < simd_width; ++j) {
		                        auto p2 = container.template at<Mask>(c2, j);
		                        f(p1, p2);
		                    }
		                }
		            }
		        }
			}

			if constexpr (vectorize) {
				for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
					AP_PREFETCH(chunks + c1 + 1);

					auto packed1 = container.template at_packed<Mask>(c1, 0);
					auto buffer1 = packed1.load_buffer();

					for (size_t j = 0; j < limit_tail; ++j) {
						auto p2 = container.template at<Mask>(c_body_end, j);

						auto buffer2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
						buffer2.force = {0,0,0};

						f(buffer1, buffer2);

						p2.force.x += buffer2.force.x.reduce_add();
						p2.force.y += buffer2.force.y.reduce_add();
						p2.force.z += buffer2.force.z.reduce_add();
					}

					packed1.force += buffer1.force;
				}
			} else {
				// body vs tail (every body chunk with tail chunk)
				for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
				    AP_PREFETCH(chunks + c1 + 1);
				    for (size_t i = 0; i < simd_width; ++i) {
				        auto p1 = container.template at<Mask>(c1, i);
				        for (size_t j = 0; j < limit_tail; ++j) {
				             auto p2 = container.template at<Mask>(c_body_end, j);
				             f(p1, p2);
				        }
				    }
				}
			}

			if constexpr (vectorize) {
				// initialize static mask elements
				alignas(64) double idx_arr[Container::chunk_size];
				for (size_t k = 0; k < Container::chunk_size; ++k) idx_arr[k] = static_cast<double>(k);
				auto lane_indices = packed::load_aligned(idx_arr);
				auto valid_tail_mask = lane_indices < static_cast<double>(limit_tail);

				const packed null = 0.0;

				// load the entire tail chunk
				auto packed_chunk = container.template at_packed<Mask>(c_body_end, 0);
				auto chunk_data = packed_chunk.load_buffer();

				// Accumulator for forces exerted by the p1 scalars onto the rest of the chunk
				pvec3 accum_force = {0, 0, 0};

				// loop over valid particles in the tail
				for (size_t i = 0; i < limit_tail; ++i) {
					auto p1 = container.template at<Mask>(c_body_end, i);
					auto buffer1 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p1);
					buffer1.force = {0, 0, 0};

					// Fresh buffer to capture forces exerted on the chunk from just this p1
					auto temp_chunk = chunk_data;
					temp_chunk.force = {0, 0, 0};

					// scalar p1 vs entire chunk
					f(buffer1, temp_chunk);

					// Dynamic mask: only interact with valid particles AHEAD of 'i' (j > i)
					// This avoids self-interaction (j == i) and double counting (j < i)
					auto current_mask = valid_tail_mask && (lane_indices > static_cast<double>(i));

					// update p1 (scalar reduction for forces ON p1 from j > i)
					p1.force += vec3 {
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
			} else {
				// tail (interact tail chunk with itself)
				for (size_t i = 0; i < limit_tail; ++i) {
  					auto p1 = container.template at<Mask>(c_body_end, i);
  					for (size_t j = i + 1; j < limit_tail; ++j) {
	   					auto p2 = container.template at<Mask>(c_body_end, j);
	   					f(p1, p2);
  					}
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










