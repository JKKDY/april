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
	concept IsBatchAtom = requires(const T& t) {
		// must have vector trait exists
		{ std::remove_cvref_t<T>::vector_trait } -> std::convertible_to<exec::internal::VectorTrait>;

		// must have a for_each_pair function
		t.template for_each_pair<ParallelPolicy::Serial, exec::internal::ExecutionMode::Hybrid>(
			universal_kernel([](auto&&, auto&&) {})
		);
	};

	template<typename T>
	concept IsBatchAtomRange = std::ranges::input_range<T> && IsBatchAtom<std::ranges::range_value_t<T>>;

	template<typename T>
	concept IsBatch = IsBatchBase<T> && (IsBatchAtom<T> || IsBatchAtomRange<T>);


	//------------
	// BCP CONCEPT
	//------------
	template<typename F>
	concept IsBCP = requires(const F& f, const vec3& v, const pvec3& pv) {
		// Must handle scalar vectors
		{ f(v) } -> std::convertible_to<vec3>;

		// Must handle packed SIMD vectors
		{ f(pv) } -> std::convertible_to<pvec3>;
	};

	struct NoBatchBCP {
		template <class T>
		constexpr T operator()(T&& v) const noexcept {
			// return std::forward<T>(v); // identity; do nothing
			return v;
		}
	};
}














