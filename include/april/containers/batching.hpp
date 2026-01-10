#pragma once
#include <cstdint>

#include "april/common.hpp"
#include "april/particle/defs.hpp"


namespace april::container {

	//-----------------
	// BATCH BASE CLASS
	//-----------------

	// batch type could also be inferred through structural matching, however setting a flag forces the user to be
	// explicit with their intent which allows for more robust concepts to enforce intent
	enum class BatchType : uint8_t {
		Symmetric,	// batch must contain symmetric work: Pairs from a single list; compute only i<j to avoid duplicate interactions.
		Asymmetric,	// batch must contain asymmetric work: Pairs from two distinct lists; compute full Cartesian product.
		Compound	// batch must contain both symmetric and asymmetric work
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

	template<BatchType batch_type, ParallelPolicy parallelize, UpdatePolicy upd, ComputePolicy cmp>
	struct BatchBase {
		static constexpr auto type = batch_type;
		static constexpr auto parallel_policy = parallelize;
		static constexpr auto update_policy = upd;
		static constexpr auto compute_policy = cmp;
		std::pair<env::ParticleType, env::ParticleType> types {};
	};

	template<BatchType type>
	using SerialBatch = BatchBase<type, ParallelPolicy::None, UpdatePolicy::Serial, ComputePolicy::Scalar>;



	namespace internal {
		// base constraints common to all batches
		template <typename T>
		concept IsBatchBase = requires(const T& b) {
			// must have static constexpr configuration flags
			{ T::type }				-> std::convertible_to<BatchType>;
			{ T::parallel_policy }	-> std::convertible_to<ParallelPolicy>;
			{ T::update_policy }	-> std::convertible_to<UpdatePolicy>;
			{ T::compute_policy }	-> std::convertible_to<ComputePolicy>;

			// must have type pair
			{ b.types } -> std::convertible_to<std::pair<env::ParticleType, env::ParticleType>>;
		};


		// a range that yields particle indices (size_t)
		template <typename T>
		concept IsIndexRange =
			std::ranges::input_range<T> &&  // must be iterable
			std::convertible_to<std::ranges::range_value_t<T>, size_t>; // iteration yields size_t


		//---------------------
		// ATOMIC UNITS OF WORK
		//---------------------
		// a single piece of work (atom). Can be either symmetric or asymmetric

		// A Symmetric Work Unit has a single range of indices (i < j logic)
		template <typename T>
		concept IsSymmetricAtom = requires(const T& t) {
			{ t.indices } -> IsIndexRange;
		};

		// An Asymmetric Work Unit has two ranges (i vs j logic)
		template <typename T>
		concept IsAsymmetricAtom = requires(const T& t) {
			{ t.indices1 } -> IsIndexRange;
			{ t.indices2 } -> IsIndexRange;
		};

		template<typename T>
		concept IsAtom = IsSymmetricAtom<T> || IsAsymmetricAtom<T>;


		//----------------------
		// RANGES OF ATOMIC WORK
		//----------------------
		// work ranges hold ranges of either symmetric or asymmetric work atoms

		// a symmetric work range contains symmetric work units
		template <typename T>
		concept IsSymmetricAtomRange =
			std::ranges::input_range<T> &&  // must be iterable
			IsSymmetricAtom<std::ranges::range_value_t<T>>;  // elements are symmetric work units

		template <typename T>
		concept HasSymmetricRange = requires(const T& b) {
			{ b.sym_chunks } -> IsSymmetricAtomRange;
		};

		// an asymmetric work range contains asymmetric work units
		template <typename T>
		concept IsAsymmetricAtomRange =
			std::ranges::input_range<T> &&  // must be iterable
			IsAsymmetricAtom<std::ranges::range_value_t<T>>;  // elements are asymmetric work units

		template <typename T>
		concept HasAsymmetricRange = requires(const T& b) {
			{ b.asym_chunks } -> IsAsymmetricAtomRange;
		};


		//-----------
		// WORK UNITS
		//-----------
		// a work unit is either an atom or a range of atoms (but cannot be both at once)
		template <typename T>
		struct SymmetricAmbiguityCheck {
			static constexpr bool value = true;
			static_assert(!(IsSymmetricAtom<T> && HasSymmetricRange<T>),
				"AMBIGUITY ERROR: Type cannot satisfy both IsSymmetricAtom and HasSymmetricRange simultaneously.");
		};

		template <typename T>
		struct AsymmetricAmbiguityCheck {
			static constexpr bool value = true;
			static_assert(!(IsAsymmetricAtom<T> && HasAsymmetricRange<T>),
				"AMBIGUITY ERROR: Type cannot satisfy both IsAsymmetricAtom and HasAsymmetricRange simultaneously.");
		};

		template <typename T>
		concept IsSymmetricWorkUnit = (IsSymmetricAtom<T> || HasSymmetricRange<T>)
			&& SymmetricAmbiguityCheck<T>::value;

		template <typename T>
		concept IsAsymmetricWorkUnit = (IsAsymmetricAtom<T> || HasAsymmetricRange<T>)
			&& AsymmetricAmbiguityCheck<T>::value;

		template<typename T> // note that T can be a symmetric work unit as well as an asymmetric one (compound unit)
		concept IsWorkUnit = IsSymmetricWorkUnit<T> || IsAsymmetricWorkUnit<T>;


		//---------------
		// COMPOUND UNITS
		//---------------
		// compound units most hold both symmetric and asymmetric work
		template <typename T>
		concept IsCompoundWorkUnit = IsSymmetricWorkUnit<T> && IsAsymmetricWorkUnit<T>;

		template <typename T>
		concept IsCompositeWorkRange =
			std::ranges::input_range<T> &&  // must be iterable
			IsCompoundWorkUnit<std::ranges::range_value_t<T>>;  // elements are compound work units

		template <typename T>
		concept HasCompositeRange = requires(const T& b) {
			{ b.chunks } -> IsCompositeWorkRange;
		};


		//---------------
		// BATCH CONCEPTS
		//---------------
		// symmetric batch: either a symmetric atom or a range of symmetric atoms
		template <typename T>
		concept IsSymmetricBatch = IsBatchBase<T> && requires {
			requires T::type == BatchType::Symmetric && IsSymmetricWorkUnit<T>;
		};

		// asymmetric batch: either an asymmetric atom or a range of asymmetric atoms
		template <typename T>
		concept IsAsymmetricBatch = IsBatchBase<T> && requires {
			requires T::type == BatchType::Asymmetric && IsAsymmetricWorkUnit<T>;
		};

		// compound batch: either a compound unit or a range of compound units
		template <typename T>
		concept IsCompoundBatch = IsBatchBase<T> && requires(const T& b) {
			requires T::type == BatchType::Compound && (IsCompoundWorkUnit<T> != HasCompositeRange<T>);
		};
	}



	//--------------
	// BATCH CONCEPT
	//--------------
	template <typename T>
	concept IsBatch =
		internal::IsAsymmetricBatch<T>	||
		internal::IsSymmetricBatch<T>	||
		internal::IsCompoundBatch<T>;



	//--------------------------
	// PREDEFINED BATCHING TYPES
	//--------------------------
	template<typename Range1, typename Range2>
	struct DirectAsymmetricBatch : SerialBatch<BatchType::Asymmetric> {
		// Holds the particle indices directly (e.g. iota_view or std::vector)
		Range1 indices1;
		Range2 indices2;
	};


	template<typename Range>
	struct DirectSymmetricBatch : SerialBatch<BatchType::Symmetric> {
		// Holds the particle indices directly (e.g. iota_view or std::vector)
		Range indices;
	};

	// Holds a list of symmetric work units (e.g. for verlet lists or grouped tasks)
	template<typename Range>
	struct ChunkedSymmetricBatch : SerialBatch<BatchType::Symmetric> {
		struct Chunk { Range indices; };

		// Satisfies HasSymmetricRange
		std::vector<Chunk> sym_chunks;
	};

	// Holds a list of asymmetric work units
	template<typename Range1, typename Range2>
	struct ChunkedAsymmetricBatch : SerialBatch<BatchType::Asymmetric> {
		struct Chunk { Range1 indices1; Range2 indices2; };

		// Satisfies HasAsymmetricRange
		std::vector<Chunk> asym_chunks;
	};

	struct TopologyBatch {
		env::ParticleID id1, id2;
		std::vector<std::pair<env::ParticleID, env::ParticleID>> pairs;
	};



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



