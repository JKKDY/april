#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"

#include "../../exec/policy.hpp"

namespace april::container::internal {

	template<typename Container>
	struct AsymmetricScalarBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::ScalarOnly> {
		explicit AsymmetricScalarBatch(Container & container) : container(container) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, typename Func>
		AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
			// Map internal ExecutionMode back to public VectorPolicy
			// constexpr VectorPolicy Pol = [] {
			// 	if constexpr (E == april::internal::ExecutionMode::Vector) {
			// 		return VectorPolicy::Vector;
			// 	} else if constexpr (E == april::internal::ExecutionMode::Scalar) {
			// 		return VectorPolicy::Scalar;
			// 	} else {
			// 		// E == ExecutionMode::Both
			// 		return VectorPolicy::Auto;
			// 	}
			// }();

			container.template for_each_particle<Mask, P, VectorPolicy::Scalar> (range1.start, range1.stop, [&](auto && p1) {
				container.template for_each_particle<Mask, P, VectorPolicy::Scalar> (range2.start, range2.stop, [&](auto && p2) {
					f(p1, p2);
				});
			});
		}

		math::Range range1;
		math::Range range2;
	private:
		Container & container;
	};

	template<typename Container>
	struct SymmetricScalarBatch : BatchBase<exec::internal::ParallelTrait::None, exec::internal::VectorTrait::ScalarOnly> {
		explicit SymmetricScalarBatch(Container & container) : container(container) {}

		template<ParticleField Mask, ParallelPolicy P, exec::internal::ExecutionMode E, typename Func>
		AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
			// TODO implement packed
			// constexpr VectorPolicy Pol = [] {
			// 	if constexpr (E == april::internal::ExecutionMode::Vector) {
			// 		return VectorPolicy::Vector;
			// 	} else if constexpr (E == april::internal::ExecutionMode::Scalar) {
			// 		return VectorPolicy::Scalar;
			// 	} else {
			// 		// E == ExecutionMode::Both
			// 		return VectorPolicy::Auto;
			// 	}
			// }();

			container.template for_each_particle<Mask, P, VectorPolicy::Scalar> (range.start, range.stop, [&](const size_t i, auto && p1) {
				container.template for_each_particle<Mask, P, VectorPolicy::Scalar> (i + 1, range.stop, [&](auto && p2) {
					f(p1, p2);
				});
			});
		}

		math::Range range;
	private:
		Container & container;
	};
}





