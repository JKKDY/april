#pragma once
#include <cstdint>
#include <utility>

#include "april/common.hpp"
#include "april/particle/defs.hpp"


namespace april::container {
	enum class BatchSymmetry : uint8_t {
		Symmetric,	// Pairs from a single list; compute only i<j to avoid duplicate interactions.
		Asymmetric  // Pairs from two distinct lists; compute full Cartesian product.
	};

	enum class ParallelPolicy : uint8_t {
		None,       // Execute immediately on the current thread (Caller owns parallelism)
	    InnerLoop,  // System spawns threads for the particle indices
	    Chunks      // System spawns threads for the chunk list
	};

	enum class UpdatePolicy : uint8_t {
		Serial, // Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring).
		Atomic  // Atomic CAS/Fetch-Add. Slower. Thread-safe for overlapping writes.
	};

	enum class ComputePolicy : uint8_t { // for future use
		Scalar,
		SIMD
	};

	template<BatchSymmetry sym, ParallelPolicy parallelize, UpdatePolicy upd>
	struct BatchBase {
		static constexpr auto symmetry = sym;
		static constexpr auto parallel_policy = parallelize;
		static constexpr auto update_policy = upd;
		std::pair<env::ParticleType, env::ParticleType> types{};
	};

	template<BatchSymmetry sym>
	using SerialBatch = BatchBase<sym, ParallelPolicy::None, UpdatePolicy::Serial>;



	// ---- Helper concepts -----

	// base constraints common to all batches
	template <typename T>
	concept IsBatchBase = requires(const T& b) {
		// must have static constexpr configuration flags
		{ T::symmetry } -> std::convertible_to<BatchSymmetry>;
		{ T::parallel_policy } -> std::convertible_to<ParallelPolicy>;
		{ T::update_policy } -> std::convertible_to<UpdatePolicy>;

		// must have type pair
		{ b.types } -> std::convertible_to<std::pair<env::ParticleType, env::ParticleType>>;
	};


	// a range that yields particle indices (size_t)
	template <typename T>
	concept IsIndexRange = std::ranges::input_range<T> &&  // must be iterable
						   std::convertible_to<std::ranges::range_value_t<T>, size_t>; // iteration yields size_t

	// checks for "indices" (Symmetric)
	template <typename T>
	concept HasSymmetricIndices = requires(const T& t) {
		{ t.indices } -> IsIndexRange;
	};

	// checks for "indices1, indices2" (Asymmetric)
	template <typename T>
	concept HasAsymmetricIndices = requires(const T& t) {
		{ t.indices1 } -> IsIndexRange;
		{ t.indices2 } -> IsIndexRange;
	};


	// ---- Main concepts -----

	// DIRECT BATCH:
	// The batch itself holds the indices.
	// if symmetric must contain a single index range
	// if asymmetric must contain two index ranges
	template <typename T>
	concept IsDirectBatch = IsBatchBase<T> && requires(const T& b) {
		requires (
		   (T::symmetry == BatchSymmetry::Symmetric  && HasSymmetricIndices<T>) ||
		   (T::symmetry == BatchSymmetry::Asymmetric && HasAsymmetricIndices<T>)
		);
	};

	// CHUNKED BATCH:
	// The batch holds a list of items (chunks), and those items hold the indices.
	template <typename T>
	concept IsChunkedBatch = IsBatchBase<T> && requires(const T& b) {
		{ b.chunks } -> std::ranges::input_range;

		requires (
			(T::symmetry == BatchSymmetry::Symmetric &&
			 HasSymmetricIndices<std::ranges::range_value_t<decltype(b.chunks)>>)
			||
			(T::symmetry == BatchSymmetry::Asymmetric &&
			 HasAsymmetricIndices<std::ranges::range_value_t<decltype(b.chunks)>>)
		);
	};

	// UNIFIED CONCEPT
	template <typename T>
	concept IsBatch = IsDirectBatch<T> || IsChunkedBatch<T>;



	// ---- boundary condition function concept -----

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



	// ---- Topology batch
	struct TopologyBatch {
		env::ParticleID id1, id2;
		std::vector<std::pair<env::ParticleID, env::ParticleID>> pairs;
	};

}
