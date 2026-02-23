#pragma once
#include <concepts>

#include "april/exec/policy.hpp"
#include "april/particle/particle_types.hpp"


namespace april::container::batching {


	//------------------------
	// CONVENIENCE DEFINITIONS
	//------------------------
	template<exec::internal::ParallelTrait P, exec::internal::VectorTrait V>
	struct BatchBase {
		static constexpr auto parallel_trait = P;
		static constexpr auto vector_trait = V;
		std::pair<ParticleType, ParticleType> types {};
	};

	struct TopologyBatch {
		ParticleID id1, id2;
		std::vector<std::pair<ParticleID, ParticleID>> pairs;
	};

	//--------------
	// BATCH CONCEPT
	//--------------
	// base constraints common to all batches
	template <typename T>
	concept IsBatchBase = requires(const T& b) {
		// must have static constexpr trait flags
		{ T::parallel_trait } -> std::convertible_to<exec::internal::ParallelTrait>;
		{ T::vector_trait }	-> std::convertible_to<exec::internal::VectorTrait>;

		// must have type pair
		{ b.types } -> std::convertible_to<std::pair<ParticleType, ParticleType>>;
	};

	template<typename T>
	concept IsBatchAtom = true ; // TODO
	// requires(const T& t) {
	// 	// callable must take in two particle views
	// 	{ t.template for_each_pair<ParticleField::all, ParallelPolicy::Serial, VectorPolicy::Scalar>(
	// 		[]<typename P0, typename P1>(P0&&, P1&&)
	// 		// requires env::IsRestrictedRef<P0> && env::IsRestrictedRef<P1>
	// 		{}
	// 	) };
	// };

	template<typename T>
	concept IsBatchAtomRange = std::ranges::input_range<T> && IsBatchAtom<std::ranges::range_value_t<T>>;

	template<typename T>
	concept IsBatch = IsBatchBase<T> && (IsBatchAtom<T> || IsBatchAtomRange<T>);


	//------------
	// BCP CONCEPT
	//------------
	template<typename F>
	concept IsBCP = true; // TODO: fix
	// requires(const F& f, const vec3& v) {
	// 	{ f(v) } -> std::convertible_to<vec3>;
	// };

	struct NoBatchBCP {
		template <class T>
		constexpr T operator()(T&& v) const noexcept {
			// return std::forward<T>(v); // identity; do nothing
			return v;
		}
	};
}














