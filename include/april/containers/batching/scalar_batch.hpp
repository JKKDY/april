#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"

namespace april::container::batching {

	template<typename Container,
		exec::VectorTrait Traits = exec::VectorTrait::ScalarPath | exec::VectorTrait::VectorPath>
	struct AsymmetricScalarBatch :
		BatchBase<exec::ParallelTrait::None, Traits> {
		explicit AsymmetricScalarBatch(Container & container) : container(container) {
			for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
		}

		template<ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Kernel>
		AP_FORCE_INLINE void for_each_pair(Kernel && f) const {
			if (range1.start == range1.stop || range2.start == range2.stop) return;

			if constexpr (static_cast<bool>(E & exec::ExecutionMode::Vector)) {
				for_each_pair_packed<P>(std::forward<Kernel>(f));
			} else {
				for_each_pair_scalar<P>(std::forward<Kernel>(f));
			}
		}

		math::Range range1;
		math::Range range2;
	private:
		Container & container;
		static constexpr size_t packed_size = packed::size();
		alignas(64) packed::value_type idx_arr[packed_size];

		// VECTORIZED EXECUTION PATH
		template<ParallelPolicy P, exec::IsKernel Kernel>
	    void for_each_pair_packed(Kernel&& f) const {
			using K = std::remove_cvref_t<Kernel>;

			// Calculate Alignment Boundaries for Range 1
			const size_t r1_rem = range1.size() % packed_size;
			const size_t tail1_start = range1.stop - r1_rem;
			const math::Range tail1 ={tail1_start,tail1_start + r1_rem};
			const math::Range body1 = {range1.start, tail1_start, packed_size};

			// Calculate Alignment Boundaries for Range 2
			const size_t r2_rem = range2.size() % packed_size;
			const size_t tail2_start = range2.stop - r2_rem;
			const math::Range tail2 ={tail2_start,tail2_start + r2_rem};
			const math::Range body2 = {range2.start, tail2_start, packed_size};

			// body1 vs body2
			for (size_t i : body1) {
				auto packed1 = container.template at_packed<K::Read, K::Write>(i);
				auto buffer1 = packed1.load_buffer();

				for (size_t j : body2) {
					auto packed2 = container.template at_packed<K::Read, K::Write>(j);
					auto buffer2 = packed2.load_buffer();
					AP_UNROLL_LOOP_N(packed_size)
					for (size_t k = 0; k < packed_size; k++) {
						auto view1 = buffer1.to_view();
						auto view2 = buffer2.to_view();
						f(view1, view2);
						buffer2.rotate_right();
					}
					buffer2.update_into(packed2);
				}
				buffer1.update_into(packed1);
			}

			// tail1 vs body2
			for (size_t i : tail1) {
				auto p1 = container.template at<K::Read, K::Write>(i);
				auto buffer1 = p1.broadcast();

				for (size_t j : body2) {
					auto packed2 = container.template at_packed<K::Read, K::Write>(j);
					auto buffer2 = packed2.load_buffer();

					auto view1 = buffer1.to_view();
					auto view2 = buffer2.to_view();
					f(view1, view2);

					buffer2.update_into(packed2);
				}

				buffer1.reduce_into(p1);
			}

			// tail2 vs body 2
			for (size_t i : tail2) {
				auto p2 = container.template at<K::Read, K::Write>(i);
				auto buffer2 = p2.broadcast();

				for (size_t j : body1) {
					auto packed1 = container.template at_packed<K::Read, K::Write>(j);
					auto buffer1 = packed1.load_buffer();

					auto view1 = buffer1.to_view();
					auto view2 = buffer2.to_view();
					f(view1, view2);

					buffer1.update_into(packed1);
				}

				buffer2.reduce_into(p2);
			}

			if (r1_rem > 0 && r2_rem > 0) {
				const auto lane_indices = packed::load_aligned(idx_arr);
				const auto mask = lane_indices < static_cast<double>(r1_rem);

				auto packed1 = container.template at_packed<K::Read, K::Write>(tail1_start);
				auto buffer1 = packed1.load_buffer();

				for (size_t j : tail2) {
					auto p2 = container.template at<K::Read, K::Write>(j);
					auto buffer2 = p2.broadcast();

					auto view1 = buffer1.to_view();
					auto view2 = buffer2.to_view();
					f(view1, view2);

					buffer2.reduce_into(p2, mask);
				}

				buffer1.update_into(packed1, mask);
			}
	    };


		// SCALAR EXECUTION PATH
		template<ParallelPolicy P, exec::IsKernel Kernel>
		void for_each_pair_scalar(Kernel&& f) const {
			using K = std::remove_cvref_t<Kernel>;
			for (size_t i = range1.start; i < range1.stop; ++i) {
				auto p1 = container.template at<K::Read, K::Write>(i);
				for (size_t j = range2.start; j < range2.stop; ++j) {
					auto p2 = container.template at<K::Read, K::Write>(j);
					f(p1, p2);
				}
			}
		}
	};





	//----------------
	// SYMMETRIC BATCH
	//----------------
	template<typename Container,
		exec::VectorTrait Traits = exec::VectorTrait::ScalarPath | exec::VectorTrait::VectorPath>
	struct SymmetricScalarBatch : BatchBase<exec::ParallelTrait::None, Traits> {
		explicit SymmetricScalarBatch(Container & container) : container(container) {
			for (size_t k = 0; k < packed_size; ++k) idx_arr[k] = static_cast<double>(k);
		}

		template<ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Kernel>
		AP_FORCE_INLINE void for_each_pair(Kernel && f) const {
			if (range.start == range.stop) return;

			if constexpr (static_cast<bool>(E & exec::ExecutionMode::Vector)) {
				for_each_pair_packed<P>(std::forward<Kernel>(f));
			} else {
				for_each_pair_scalar<P>(std::forward<Kernel>(f));
			}
		}

		math::Range range;
	private:
		Container & container;
		static constexpr size_t packed_size = packed::size();
		alignas(64) packed::value_type idx_arr[packed_size];


		// VECTORIZED EXECUTION PATH
		template<ParallelPolicy P, exec::IsKernel Kernel>
	    void for_each_pair_packed(Kernel&& f) const {
	        using K = std::remove_cvref_t<Kernel>;

	        // Calculate Alignment Boundaries
	        const size_t rem = range.size() % packed_size;
	        const size_t tail_start = range.stop - rem;
	        const math::Range tail = {tail_start, range.stop};
	        const math::Range body = {range.start, tail_start, packed_size};

	        // body vs body (Upper Triangle)
	        for (size_t i : body) {
	            auto packed1 = container.template at_packed<K::Read, K::Write>(i);
	            auto buffer1 = packed1.load_buffer();

	            // self-interaction within block i
	            {
	                auto buffer2 = packed1.load_buffer();
	                AP_UNROLL_LOOP()
	                for (size_t k = 0; k < packed_size / 2 - 1; k++) {
	                    buffer2.rotate_right();
	                    auto view1 = buffer1.to_view();
	                    auto view2 = buffer2.to_view();
	                    f(view1, view2);
	                }

	                buffer2.template rotate_left<packed_size / 2 - 1>();
	                buffer1.accumulate(buffer2);
	                buffer2.template rotate_right<packed_size / 2>();

	                auto view1 = buffer1.to_view();
	                auto view2 = buffer2.to_view();
	                f(view1, view2);
	            }

	            // block i vs block j (where j > i)
	            for (size_t j = i + packed_size; j < tail_start; j += packed_size) {
	                auto packed2 = container.template at_packed<K::Read, K::Write>(j);
	                auto buffer2 = packed2.load_buffer();
	                AP_UNROLL_LOOP_N(packed_size)
	                for (size_t k = 0; k < packed_size; k++) {
	                    auto view1 = buffer1.to_view();
	                    auto view2 = buffer2.to_view();
	                    f(view1, view2);
	                    buffer2.rotate_right();
	                }
	                buffer2.update_into(packed2);
	            }
	            buffer1.update_into(packed1);
	        }

	        // tail vs body
	        for (size_t i : tail) {
	            auto p1 = container.template at<K::Read, K::Write>(i);
	            auto buffer1 = p1.broadcast();

	            for (size_t j : body) {
	                auto packed2 = container.template at_packed<K::Read, K::Write>(j);
	                auto buffer2 = packed2.load_buffer();

	                auto view1 = buffer1.to_view();
	                auto view2 = buffer2.to_view();
	                f(view1, view2);

	                buffer2.update_into(packed2);
	            }
	            buffer1.reduce_into(p1);
	        }

	        // tail vs tail (Masked Upper Triangle)
	        if (rem > 0) {
	            const auto absolute_lane_indices = packed::load_aligned(idx_arr) + static_cast<double>(tail_start);
	            const auto valid_tail_mask = absolute_lane_indices < static_cast<double>(range.stop);

	            auto packed_chunk = container.template at_packed<K::Read, K::Write>(tail_start);
	            auto chunk_data = packed_chunk.load_buffer();

	            for (size_t i : tail) {
	                auto p1 = container.template at<K::Read, K::Write>(i);
	                auto buffer1 = p1.broadcast();

	                // Temp buffer to isolate forces applied from p1
	                auto temp_chunk = packed_chunk.load_buffer();

	                auto view1 = buffer1.to_view();
	                auto view2 = temp_chunk.to_view();
	                f(view1, view2);

	                auto current_mask = valid_tail_mask && (absolute_lane_indices > static_cast<double>(i));
	                buffer1.reduce_into(p1, current_mask);
	                chunk_data.accumulate(temp_chunk, current_mask);
	            }
	            chunk_data.update_into(packed_chunk);
	        }
	    }


		// SCALAR EXECUTION PATH
		template<ParallelPolicy P, exec::IsKernel Kernel>
		void for_each_pair_scalar(Kernel&& f) const {
			using K = std::remove_cvref_t<Kernel>;
			for (size_t i = range.start; i < range.stop; ++i) {
				auto p1 = container.template at<K::Read, K::Write>(i);
				for (size_t j = i + 1; j < range.stop; ++j) {
					auto p2 = container.template at<K::Read, K::Write>(j);
					f(p1, p2);
				}
			}
		}
	};
}
















