#pragma once
#include <cstdint>
#include <concepts>

#include "april/base/types.hpp"
#include "april/particle/defs.hpp"


namespace april::container {

	//---------------
	// BATCH POLICIES
	//---------------
	enum class ParallelTrait : uint8_t {
		None,       // Execute immediately on the current thread (Caller owns parallelism)
	    Inner,		// System spawns threads for executing a single batch
	};

	enum class UpdatePolicy : uint8_t {
		Serial,		// Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring).
		Atomic,		// Atomic CAS/Fetch-Add. Slower. Thread-safe for overlapping writes.
	};

	enum class ComputeTrait : uint8_t {
		Scalar,
		Vector,
		Hybrid
	};


	//------------------------
	// CONVENIENCE DEFINITIONS
	//------------------------
	template<ParallelTrait parallelize, UpdatePolicy upd, ComputeTrait cmp>
	struct BatchBase {
		static constexpr auto parallel_policy = parallelize;
		static constexpr auto update_policy = upd;
		static constexpr auto compute_policy = cmp;
		std::pair<ParticleType, ParticleType> types {};
	};

	using SerialBatch = BatchBase<ParallelTrait::None, UpdatePolicy::Serial, ComputeTrait::Scalar>;

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
		// must have static constexpr configuration flags
		{ T::parallel_policy }	-> std::convertible_to<ParallelTrait>;
		{ T::update_policy }	-> std::convertible_to<UpdatePolicy>;
		{ T::compute_policy }	-> std::convertible_to<ComputeTrait>;

		// must have type pair
		{ b.types } -> std::convertible_to<std::pair<ParticleType, ParticleType>>;
	};

	template<typename T>
	concept IsBatchAtom = requires(const T& t) {
		// callable must take in two particle views
		{ t.template for_each_pair<ParticleField::all>(
			[]<typename P0, typename P1>(P0&&, P1&&)
			// requires env::IsRestrictedRef<P0> && env::IsRestrictedRef<P1>
			{}
		) };
	};

	template<typename T>
	concept IsBatchAtomRange = std::ranges::input_range<T> && IsBatchAtom<std::ranges::range_value_t<T>>;

	template<typename T>
	concept IsBatch = IsBatchBase<T> && (IsBatchAtom<T> || IsBatchAtomRange<T>);


	//------------
	// BCP CONCEPT
	//------------
	template<typename F>
	concept IsBCP = requires(const F& f, const vec3& v) {
		{ f(v) } -> std::convertible_to<vec3>;
	};

	struct NoBatchBCP {
		template <class T>
		constexpr T operator()(T&& v) const noexcept {
			return std::forward<T>(v); // identity; do nothing
		}
	};
}





