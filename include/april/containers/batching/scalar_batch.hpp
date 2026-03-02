#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"

namespace april::container::batching {

	template<typename Container>
	struct AsymmetricScalarBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::ScalarOnly> {
		explicit AsymmetricScalarBatch(Container & container) : container(container) {}

		template<ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Kernel>
		AP_FORCE_INLINE
		void for_each_pair (Kernel && f) const {
			using K = std::remove_cvref_t<Kernel>;
			container.template for_each_particle<P, VectorPolicy::Scalar> (range1.start, range1.stop, scalar_kernel<K::Read, K::Write>(
				[&]<bool is_packed1>(auto && p1) {
					container.template for_each_particle<P,  VectorPolicy::Scalar> (range2.start, range2.stop, scalar_kernel<K::Read, K::Write>(
						[&]<bool is_packed2>(auto && p2) {
							f(p1, p2);
						})
					);
				})
			);
		}

		math::Range range1;
		math::Range range2;
	private:
		Container & container;
	};

	template<typename Container>
	struct SymmetricScalarBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::ScalarOnly> {
		explicit SymmetricScalarBatch(Container & container) : container(container) {}

		template<ParallelPolicy P, exec::internal::ExecutionMode E, exec::IsKernel Kernel>
		AP_FORCE_INLINE
		void for_each_pair (Kernel && f) const {
			using K = std::remove_cvref_t<Kernel>;
			container.template for_each_particle<P, VectorPolicy::Scalar> (range.start, range.stop, scalar_kernel<K::Read, K::Write>(
				[&]<bool is_packed1>(const size_t i, auto && p1) {
					container.template for_each_particle<P, VectorPolicy::Scalar> (i + 1, range.stop, scalar_kernel<K::Read, K::Write>(
						[&]<bool is_packed2>(auto && p2) {
							f(p1, p2);
						})
					);
				})
			);
		}

		math::Range range;
	private:
		Container & container;
	};
}














