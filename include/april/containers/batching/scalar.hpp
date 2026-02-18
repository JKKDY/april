#pragma once

#include "april/base/macros.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"
#include "april/particle/fields.hpp"
#include "april/base/policy.hpp"

namespace april::container::internal {

	template<typename Container>
	struct AsymmetricScalarBatch : SerialBatch {
		explicit AsymmetricScalarBatch(Container & container) : container(container) {}

		template<ParticleField Mask, ParallelPolicy P, VectorPolicy V, typename Func>
		AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
			container.template for_each_particle<Mask, P, V> (range1.start, range1.stop, [&](auto && p1) {
				container.template for_each_particle<Mask, P, V> (range2.start, range2.stop, [&](auto && p2) {
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
	struct SymmetricScalarBatch : SerialBatch {
		explicit SymmetricScalarBatch(Container & container) : container(container) {}

		template<ParticleField Mask, ParallelPolicy P, VectorPolicy V, typename Func>
		AP_FORCE_INLINE
		void for_each_pair (Func && f) const {
			container.template for_each_particle<Mask, P, V> (range.start, range.stop, [&](const size_t i, auto && p1) {
				container.template for_each_particle<Mask, P, V> (i + 1, range.stop, [&](auto && p2) {
					f(p1, p2);
				});
			});
		}

		math::Range range;
	private:
		Container & container;
	};
}





