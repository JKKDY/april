#pragma once
#include <cstdint>
#include <concepts>

#include "april/base/policy.hpp"
#include "april/particle/defs.hpp"


namespace april::container {

	//---------------
	// BATCH POLICIES
	//---------------
	// enum class ParallelTrait : uint8_t {
	// 	None,       // Execute immediately on the current thread (Caller owns parallelism)
	//     Inner,		// System spawns threads for executing a single batch
	// };

	// enum class UpdatePolicy : uint8_t {
	// 	Serial,		// Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring).
	// 	Atomic,		// Atomic CAS/Fetch-Add. Slower. Thread-safe for overlapping writes.
	// };

	// enum class ComputeTrait : uint8_t {
	// 	Scalar,
	// 	Vector,
	// 	Hybrid
	// };


	//------------------------
	// CONVENIENCE DEFINITIONS
	//------------------------
	template<april::internal::ParallelTrait P, april::internal::VectorTrait V>
	struct BatchBase {
		static constexpr auto parallel_trait = P;
		static constexpr auto vector_trait = V;
		std::pair<ParticleType, ParticleType> types {};
	};

	using SerialBatch = BatchBase<april::internal::ParallelTrait::None, april::internal::VectorTrait::ScalarOnly>;

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
		{ T::parallel_trait } -> std::convertible_to< april::internal::ParallelTrait>;
		{ T::vector_trait }	-> std::convertible_to< april::internal::VectorTrait>;

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





