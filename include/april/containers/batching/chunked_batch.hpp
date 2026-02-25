#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"


namespace april::container::batching {

   	template<typename Container, typename ChunkPtr>
	struct AsymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
		explicit AsymmetricChunkedBatch(Container & container, ChunkPtr * chunks): container(container), chunks(chunks) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
		AP_FORCE_INLINE
   		void for_each_pair (Func && f) const {
			// skip empty range
		    if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

			if constexpr (static_cast<bool>(E & exec::internal::ExecutionMode::Vector)) {
				for_each_pair_packed<Mask, P>(f);
			} else {
				for_each_pair_scalar<Mask, P>(f);
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

   		static constexpr size_t chunk_size = Container::chunk_size;
   		static constexpr size_t packed_size = packed::size();
   		static constexpr  size_t iter_chunks = chunk_size/packed_size;


   		template<ParticleField Mask, ParallelPolicy P>
		void for_each_pair_packed(auto && f) const {
   			// handle full chunks masquerading as 0-length tails
   			const size_t tail1 = (range1_tail == 0) ? chunk_size : range1_tail;
   			const size_t tail2 = (range2_tail == 0) ? chunk_size : range2_tail;

   			// calculate exact index limits
   			const size_t full_tail_end1 = (tail1 / packed_size) * packed_size;
   			const size_t full_tail_end2 = (tail2 / packed_size) * packed_size;

   			// Define the ranges
   			// range of full chunks (tail chunks peeled off).Iterates over chunks
   			const math::Range full_chunks1 = {range1_chunks.start, range1_chunks.stop - 1};
   			const math::Range full_chunks2 = {range2_chunks.start, range2_chunks.stop - 1};

   			// range of full simd widths inside tail chunks. Iterates in packed_size steps
   			const math::Range full_tail1 = {0, full_tail_end1, packed_size};
   			const math::Range full_tail2 = {0, full_tail_end2, packed_size};

   			// number of particles at the end of each tail-tail. Iterates over individual particles
   			const math::Range partial_tail1 = {full_tail_end1, tail1};
   			const math::Range partial_tail2 = {full_tail_end2, tail2};


   			//-------------
   			// BODY VS BODY
   			//-------------
   			// helper for full n x n kernel
   			auto interact = [&](auto & buffer1, auto & buffer2) {
   				AP_UNROLL_LOOP_N(packed_size)
				   for (size_t k = 0; k < packed_size; k++) {
				   	f(buffer1, buffer2);
				   	buffer2.rotate_right();
				   }
   			};

   			// lambda to sweep a hoisted buffer1 across all full SIMD blocks in range 2
   			auto sweep_full_range2 = [&](auto & buffer1) {
   				for (size_t c2 : full_chunks2) {
   					AP_PREFETCH(chunks + c2 + 1);
   					AP_UNROLL_LOOP_N(iter_chunks)
					   for (size_t j = 0; j < chunk_size; j += packed_size) {
					   	auto packed2 = container.template at_packed<Mask>(c2, j);
					   	auto buffer2 = packed2.load_buffer();
					   	interact(buffer1, buffer2);
					   	packed2.force = buffer2.force;
					   }
   				}
   				for (size_t t2 : full_tail2) {
   					auto packed2 = container.template at_packed<Mask>(full_chunks2.stop, t2);
   					auto buffer2 = packed2.load_buffer();
   					interact(buffer1, buffer2);
   					packed2.force = buffer2.force;
   				}
   			};

		    // full Chunks 1 vs Full Range 2 ()
		    for (size_t c1 : full_chunks1) {
		        AP_PREFETCH(chunks + c1 + 1);
		        AP_UNROLL_LOOP_N(iter_chunks)
		        for (size_t i = 0; i < chunk_size; i += packed_size) {
		            auto packed1 = container.template at_packed<Mask>(c1, i);
		            auto buffer1 = packed1.load_buffer();

		            sweep_full_range2(buffer1);
		            packed1.force = buffer1.force;
		        }
		    }
		    // Full Tail 1 vs Full Range 2 (chunks + full tail)
		    for (size_t t1 : full_tail1) {
		        auto packed1 = container.template at_packed<Mask>(full_chunks1.stop, t1);
		        auto buffer1 = packed1.load_buffer();

		        sweep_full_range2(buffer1);
		        packed1.force = buffer1.force;
		    }


   			//--------------
   			// TAILS VS BODY
   			//--------------
   			//lambda to sweep across all full SIMD blocks in a given range
   			auto sweep_packed_blocks = [&](const math::Range& full_chunks, const math::Range& full_tail, size_t tail_chunk_idx, auto&& kernel) {
   				for (size_t c : full_chunks) {
   					AP_PREFETCH(chunks + c + 1);
   					AP_UNROLL_LOOP_N(chunk_size/packed_size)
					   for (size_t j = 0; j < chunk_size; j += packed_size) {
					   	kernel(c, j);
					   }
   				}
   				for (size_t t : full_tail) {
   					kernel(tail_chunk_idx, t);
   				}
		    };

   			// Partial Tail 2 vs Full Range 1 (chunks + full tail)
   			for (size_t i : partial_tail2) {
   				auto p2 = container.template at<Mask>(full_chunks2.stop, i);
   				auto buffer2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
   				buffer2.force = {0,0,0};

   				sweep_packed_blocks(full_chunks1, full_tail1, full_chunks1.stop, [&](size_t c1, size_t j) {
				   auto packed1 = container.template at_packed<Mask>(c1, j);
				   auto buffer1 = packed1.load_buffer();
				   buffer1.force = {0,0,0};

				   f(buffer1, buffer2);

				   packed1.force += buffer1.force;
			   });

   				p2.force.x += buffer2.force.x.reduce_add();
   				p2.force.y += buffer2.force.y.reduce_add();
   				p2.force.z += buffer2.force.z.reduce_add();
   			}

   			// Partial Tail 1 vs Full Range 2 (chunks + full tail)
   			for (size_t i : partial_tail1) {
   				auto p1 = container.template at<Mask>(full_chunks1.stop, i);
   				auto buffer1 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p1);
   				buffer1.force = {0,0,0};

   				sweep_packed_blocks(full_chunks2, full_tail2, full_chunks2.stop, [&](size_t c2, size_t j) {
				   auto packed2 = container.template at_packed<Mask>(c2, j);
				   auto buffer2 = packed2.load_buffer();
				   buffer2.force = {0,0,0};

				   f(buffer1, buffer2);

				   packed2.force += buffer2.force;
			   });

   				p1.force.x += buffer1.force.x.reduce_add();
   				p1.force.y += buffer1.force.y.reduce_add();
   				p1.force.z += buffer1.force.z.reduce_add();
   			}

   			//-------------
   			// TAIL VS TAIL
   			//-------------
   			// partial tail vs partial tail
   			if (partial_tail1.start != partial_tail1.stop && partial_tail2.start != partial_tail2.stop) {

   				// initialize mask for the valid lanes in partial_tail1
   				alignas(64) double idx_arr[packed_size];
   				for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
   				const auto lane_indices = packed::load_aligned(idx_arr);

   				const double valid_lanes = static_cast<double>(partial_tail1.stop - partial_tail1.start);
   				const auto mask = lane_indices < valid_lanes;
   				const packed null = 0.0;

   				// load the single SIMD block containing partial_tail1
   				auto packed1 = container.template at_packed<Mask>(full_chunks1.stop, partial_tail1.start);
   				auto buffer1 = packed1.load_buffer();
   				buffer1.force = {0, 0, 0};

   				// loop over the individual particles in partial_tail2
   				for (size_t i = partial_tail2.start; i < partial_tail2.stop; ++i) {
   					auto p2 = container.template at<Mask>(full_chunks2.stop, i);
   					auto buffer2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
   					buffer2.force = {0, 0, 0};

   					f(buffer1, buffer2);

   					p2.force += vec3 {
   						select(mask, buffer2.force.x, null).reduce_add(),
   						select(mask, buffer2.force.y, null).reduce_add(),
   						select(mask, buffer2.force.z, null).reduce_add()
   					};
   				}

   				// write back valid lanes to memory for packed1
   				packed1.force += pvec3 {
   					select(mask, buffer1.force.x, null),
   					select(mask, buffer1.force.y, null),
   					select(mask, buffer1.force.z, null)
   				};
   			}
		}

   		template<ParticleField Mask, ParallelPolicy P>
   		void for_each_pair_packed2 (auto && f) const {
   						// skip empty range
		    if (range1_chunks.start == range1_chunks.stop || range2_chunks.start == range2_chunks.stop) return;

		    constexpr size_t simd_width = Container::chunk_size;

			static_assert(packed::size() == Container::chunk_size);

			// peel of last chunk (i.e. the tail)
		    const size_t c1_body_end = range1_chunks.stop - 1;
		    const size_t c2_body_end = range2_chunks.stop - 1;

		    // if tail is 0, it means the chunk is actually full (stride)
		    const size_t limit1_tail = (range1_tail == 0) ? simd_width : range1_tail;
		    const size_t limit2_tail = (range2_tail == 0) ? simd_width : range2_tail;

		    // body vs body with hardcoded stride x stride loop
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


			// body 1 vs tail 2 (iterate full chunks of R1 against the single partial chunk of R2)
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

			// body 2 vs tail 1 (iterate full chunks of R2 against the single partial chunk of R1
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
   		}
	};




	template<typename Container, typename ChunkPtr>
	struct SymmetricChunkedBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::Mixed> {
		explicit SymmetricChunkedBatch(Container & container, ChunkPtr * chunks) : container(container), chunks(chunks) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Func>
	    AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
	        if (range_chunks.start == range_chunks.stop) return;

			if constexpr (static_cast<bool>(E & exec::internal::ExecutionMode::Vector)) {
				for_each_pair_scalar<Mask, P>(f);
			} else {
				for_each_pair_scalar<Mask, P>(f);
			}
		}

		// Range represents chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
		math::Range  range_chunks;  // Chunk Indices [start, end)
		size_t range_tail{};  // Number of valid items in the last chunk of range1 (0 = Full)
	private:
		Container & container;
		ChunkPtr * AP_RESTRICT const chunks = container.ptr_chunks;

		static constexpr size_t chunk_size = Container::chunk_size;
		static constexpr size_t packed_size = packed::size();

		template<ParticleField Mask, ParallelPolicy P>
		void for_each_pair_scalar (auto && f) const {
			const size_t c_body_end = range_chunks.stop - 1;
			const size_t limit_tail = (range_tail == 0) ? packed_size : range_tail;

			for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
				AP_PREFETCH(chunks + c1 + 1);

				// chunk self interaction
				for (size_t i = 0; i < chunk_size; ++i) {
					auto p1 = container.template at<Mask>(c1, i);
					for (size_t j = i + 1; j < chunk_size; ++j) {
						auto p2 = container.template at<Mask>(c1, j);
						f(p1, p2);
					}
				}

				// interaction with all other chunks
				for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
					AP_PREFETCH(chunks + c2 + 1);
					for (size_t i = 0; i < chunk_size; ++i) {
						auto p1 = container.template at<Mask>(c1, i);
						for (size_t j = 0; j < chunk_size; ++j) {
							auto p2 = container.template at<Mask>(c2, j);
							f(p1, p2);
						}
					}
				}
			}

			// body vs tail (every body chunk with tail chunk)
			for (size_t c1 = range_chunks.start; c1 < c_body_end; ++c1) {
				AP_PREFETCH(chunks + c1 + 1);
				for (size_t i = 0; i < chunk_size; ++i) {
					auto p1 = container.template at<Mask>(c1, i);
					for (size_t j = 0; j < limit_tail; ++j) {
						auto p2 = container.template at<Mask>(c_body_end, j);
						f(p1, p2);
					}
				}
			}

			// tail (interact tail chunk with itself)
			for (size_t i = 0; i < limit_tail; ++i) {
				auto p1 = container.template at<Mask>(c_body_end, i);
				for (size_t j = i + 1; j < limit_tail; ++j) {
					auto p2 = container.template at<Mask>(c_body_end, j);
					f(p1, p2);
				}
			}
		}

		template<ParticleField Mask, ParallelPolicy P>
		void for_each_pair_packed (auto && f) const {
			const size_t c_body_end = range_chunks.stop - 1;
			const size_t limit_tail = (range_tail == 0) ? packed_size : range_tail;

			// body (iterate c1 up to the last full chunk)
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
					for (size_t i = 0; i < packed_size/2 - 1; i++) {
						buffer2.rotate_right();
						f(buffer1, buffer2);
					}

					packed1.force.x += buffer2.force.x.template rotate_left<packed_size/2 - 1>();
					packed1.force.y += buffer2.force.y.template rotate_left<packed_size/2 - 1>();
					packed1.force.z += buffer2.force.z.template rotate_left<packed_size/2 - 1>();

					buffer2.rotate_right();
					f(buffer1, buffer2);
				}

				auto packed1 = container.template at_packed<Mask>(c1, 0);
				auto buffer1 = packed1.load_buffer();

				for (size_t c2 = c1 + 1; c2 < c_body_end; ++c2) {
					AP_PREFETCH(chunks + c2 + 1);
					auto packed2 = container.template at_packed<Mask>(c2, 0);
					auto buffer2 = packed2.load_buffer();

					AP_UNROLL_LOOP_N(packed_size)
					for (size_t k = 0; k < packed_size; k++) {
						f(buffer1, buffer2);
						buffer2.template rotate_right<1>();
					}
					packed2.force = buffer2.force;
				}
				packed1.force = buffer1.force;
			}

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

			// const size_t tail_len = (range_tail == 0) ? chunk_size : range_tail;
	  //       const size_t full_tail_end = (tail_len / packed_size) * packed_size;
	  //
	  //       const math::Range full_chunks = {range_chunks.start, range_chunks.stop - 1};
	  //       const math::Range full_tail = {0, full_tail_end, packed_size};
	  //       const math::Range partial_tail = {full_tail_end, tail_len};
	  //
	  //       const size_t c_tail = full_chunks.stop;
	  //
	  //       // Lambda: Computes block-to-block interaction
	  //       auto interact_blocks = [&](auto& b1, auto& b2) {
	  //           AP_UNROLL_LOOP_N(packed_size)
	  //           for (size_t k = 0; k < packed_size; k++) {
	  //               f(b1, b2);
	  //               b2.rotate_right();
	  //           }
	  //       };
	  //
	  //       // Lambda: Self-interacts a block and returns the live register
	  //       auto interact_self = [&](size_t c, size_t i) {
	  //           auto packed1 = container.template at_packed<Mask>(c, i);
	  //           auto b1 = packed1.load_buffer(); // Contains existing memory forces
	  //
	  //           auto b2 = packed1.load_buffer();
	  //           b2.force = {0,0,0}; // Tracks only the rotated deltas
	  //
	  //           AP_UNROLL_LOOP()
	  //           for (size_t r = 0; r < packed_size/2 - 1; r++) {
	  //               b2.rotate_right();
	  //               f(b1, b2);
	  //           }
	  //
	  //           b1.force.x += b2.force.x.template rotate_left<packed_size/2 - 1>();
	  //           b1.force.y += b2.force.y.template rotate_left<packed_size/2 - 1>();
	  //           b1.force.z += b2.force.z.template rotate_left<packed_size/2 - 1>();
	  //
	  //           b2.rotate_right();
	  //           f(b1, b2);
	  //
	  //           return b1; // Do NOT write to memory here! Return the accumulator.
	  //       };
	  //
	  //       // 1. Process Full Chunks
	  //       for (size_t c1 : full_chunks) {
	  //           AP_PREFETCH(chunks + c1 + 1);
	  //
	  //           for (size_t i = 0; i < chunk_size; i += packed_size) {
	  //
	  //               // Get the accumulator with the self-interaction already applied
	  //               auto b1 = interact_self(c1, i);
	  //
	  //               // Interact with j > i inside the SAME chunk
	  //               for (size_t j = i + packed_size; j < chunk_size; j += packed_size) {
	  //                   auto packed2 = container.template at_packed<Mask>(c1, j);
	  //                   auto b2 = packed2.load_buffer();
	  //                   interact_blocks(b1, b2);
	  //                   packed2.force = b2.force;
	  //               }
	  //
	  //               // Interact with c2 > c1 chunks
	  //               for (size_t c2 = c1 + 1; c2 < full_chunks.stop; ++c2) {
	  //                   for (size_t j = 0; j < chunk_size; j += packed_size) {
	  //                       auto packed2 = container.template at_packed<Mask>(c2, j);
	  //                       auto b2 = packed2.load_buffer();
	  //                       interact_blocks(b1, b2);
	  //                       packed2.force = b2.force;
	  //                   }
	  //               }
	  //
	  //               // Interact with full tail blocks
	  //               for (size_t t : full_tail) {
	  //                   auto packed_t = container.template at_packed<Mask>(c_tail, t);
	  //                   auto b2 = packed_t.load_buffer();
	  //                   interact_blocks(b1, b2);
	  //                   packed_t.force = b2.force;
	  //               }
	  //
	  //               // Interact with partial tail particles
	  //               for (size_t t : partial_tail) {
	  //                   auto p2 = container.template at<Mask>(c_tail, t);
	  //                   auto b2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
	  //                   b2.force = {0,0,0};
	  //                   f(b1, b2);
	  //                   p2.force.x += b2.force.x.reduce_add();
	  //                   p2.force.y += b2.force.y.reduce_add();
	  //                   p2.force.z += b2.force.z.reduce_add();
	  //               }
	  //
	  //               // Finally write b1 back to memory ONCE
	  //               auto packed1 = container.template at_packed<Mask>(c1, i);
	  //               packed1.force = b1.force;
	  //           }
	  //       }
	  //
	  //       // 2. Process Tail Chunk Self Interactions
	  //       for (size_t i : full_tail) {
	  //           auto b1 = interact_self(c_tail, i);
	  //
	  //           for (size_t j = i + packed_size; j < full_tail_end; j += packed_size) {
	  //               auto packed2 = container.template at_packed<Mask>(c_tail, j);
	  //               auto b2 = packed2.load_buffer();
	  //               interact_blocks(b1, b2);
	  //               packed2.force = b2.force;
	  //           }
	  //
	  //           for (size_t j : partial_tail) {
	  //               auto p2 = container.template at<Mask>(c_tail, j);
	  //               auto b2 = particle::internal::PackedParticleBuffer<Mask>::broadcast(p2);
	  //               b2.force = {0,0,0};
	  //               f(b1, b2);
	  //               p2.force.x += b2.force.x.reduce_add();
	  //               p2.force.y += b2.force.y.reduce_add();
	  //               p2.force.z += b2.force.z.reduce_add();
	  //           }
	  //
	  //           auto packed1 = container.template at_packed<Mask>(c_tail, i);
	  //           packed1.force = b1.force;
	  //       }
	  //
	  //       // 3. Partial tail vs Partial tail
	  //       // Scalar fallback: It's trivial (at most 21 interactions) and avoids SIMD permutation hell
	  //       for (size_t i = partial_tail.start; i < partial_tail.stop; ++i) {
	  //           auto p1 = container.template at<Mask>(c_tail, i);
	  //           for (size_t j = i + 1; j < partial_tail.stop; ++j) {
	  //               auto p2 = container.template at<Mask>(c_tail, j);
	  //               f(p1, p2);
	  //           }
	  //       }
		}
	};

}

