#pragma once
#include <cstdint>
#include <concepts>

#include "april/macros.hpp"
#include "april/common.hpp"
#include "april/particle/defs.hpp"
#include "april/particle/access.hpp"


namespace april::batching {

	//---------------
	// BATCH POLICIES
	//---------------
	enum class ParallelPolicy : uint8_t {
		None,       // Execute immediately on the current thread (Caller owns parallelism)
	    Inner,		// System spawns threads for executing a single batch
	};

	enum class UpdatePolicy : uint8_t {
		Serial,		// Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring).
		Serial_N3,	// Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring). Uses newton 3 to update opposite particle as well
		Atomic,		// Atomic CAS/Fetch-Add. Slower. Thread-safe for overlapping writes.
		Atomic_N3,	// Standard '+='. Fastest. Assumes thread-safety (Serial or Coloring). Uses newton 3 to update opposite particle as well
	};

	enum class ComputePolicy : uint8_t { // for future use
		Scalar,
		Vector
	};


	//------------------------
	// CONVENIENCE DEFINITIONS
	//------------------------
	template<ParallelPolicy parallelize, UpdatePolicy upd, ComputePolicy cmp>
	struct BatchBase {
		static constexpr auto parallel_policy = parallelize;
		static constexpr auto update_policy = upd;
		static constexpr auto compute_policy = cmp;
		std::pair<env::ParticleType, env::ParticleType> types {};
	};

	using SerialBatch = BatchBase<ParallelPolicy::None, UpdatePolicy::Serial, ComputePolicy::Scalar>;

	struct TopologyBatch {
		env::ParticleID id1, id2;
		std::vector<std::pair<env::ParticleID, env::ParticleID>> pairs;
	};

	//--------------
	// BATCH CONCEPT
	//--------------
	// base constraints common to all batches
	template <typename T>
	concept IsBatchBase = requires(const T& b) {
		// must have static constexpr configuration flags
		{ T::parallel_policy }	-> std::convertible_to<ParallelPolicy>;
		{ T::update_policy }	-> std::convertible_to<UpdatePolicy>;
		{ T::compute_policy }	-> std::convertible_to<ComputePolicy>;

		// must have type pair
		{ b.types } -> std::convertible_to<std::pair<env::ParticleType, env::ParticleType>>;
	};

	template<typename T>
	concept IsBatchAtom = requires(const T& t) {
		// callable must take in two particle views
		{ t.template for_each_pair<+env::Field::all>(
			[]<typename P0, typename P1>(P0&&, P1&&)
			requires env::IsRestrictedRef<P0> && env::IsRestrictedRef<P1>
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


	//--------------------
	// PRE DEFINED BATCHES
	//--------------------

	struct Range {
		Range() = default;

		// From a pair of integers
		template<std::integral I>
		Range(const std::pair<I, I>& r)
			: Range(static_cast<size_t>(r.first), static_cast<size_t>(r.second)) {}

		// explicit start/end
		Range(const size_t start, const size_t end)
			: start(start), end(end), size(end-start), empty(size==0) {}

		size_t start{}, end{};
		size_t size{};
		bool empty{};
	};

	template<typename Container>
	struct AsymmetricScalarBatch : SerialBatch {
		explicit AsymmetricScalarBatch(Container & container) : container(container) {}

		template<env::FieldMask Mask, typename Func>
		AP_FORCE_INLINE void for_each_pair (Func && f) const {
			for (size_t i = range1.start; i < range1.end; ++i) {
				auto p1 = container.template restricted_at<Mask>(i);

				for (size_t j = range2.start; j < range2.end; ++j) {
					auto p2 = container.template restricted_at<Mask>(j);
					f(p1, p2);
				}
			}
		}

		Range range1;
		Range range2;
	private:
		Container & container;
	};

	template<typename Container>
	struct SymmetricScalarBatch : SerialBatch {
		explicit SymmetricScalarBatch(Container & container) : container(container) {}

		template<env::FieldMask Mask, typename Func>
		AP_FORCE_INLINE void for_each_pair (Func && f) const {
			for (size_t i = range.start; i < range.end; ++i) {
				auto p1 = container.template restricted_at<Mask>(i);

				for (size_t j = i + 1; j < range.end; ++j) {
					auto p2 = container.template restricted_at<Mask>(j);
					f(p1, p2);
				}
			}
		}

		Range range;
	private:
		Container & container;
	};


	template<typename Container>
	requires requires(Container& c) {
		{ Container::chunk_size } -> std::convertible_to<size_t>;
		{ c.template restricted_at<+env::Field::position>(size_t{0}, size_t{0}) };
	}
	struct AsymmetricChunkedBatch : SerialBatch {
		explicit AsymmetricChunkedBatch(Container & container) : container(container) {}

		template<env::FieldMask Mask, typename Func>
		AP_FORCE_INLINE void for_each_pair (Func && f) const {
			// loop over chunks
			for (size_t c1 = range1.start; c1 < range1.end; ++c1) {
				for (size_t c2 = range2.start; c2 < range2.end; ++c2) {

					// loop inside chunks
					for (size_t i = 0; i < Container::chunk_size; ++i) {
						auto p1 = container.template restricted_at<Mask>(c1, i);
						for (size_t j = 0; j < Container::chunk_size; ++j) {
							auto p2 = container.template restricted_at<Mask>(c2, j);
							f(p1, p2);
						}
					}
				}
			}
		}


		// Ranges represent chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
		Range range1;
		Range range2;
	private:
		Container & container;
	};


	template<typename Container>
	requires requires(Container& c) {
		// Same constraints
		{ Container::chunk_size } -> std::convertible_to<size_t>;
		{ c.template restricted_at<env::Field::position>(size_t{0}, size_t{0}) };
	}
	struct SymmetricChunkedBatch : SerialBatch {
		explicit SymmetricChunkedBatch(Container & container) : container(container) {}

		template<env::FieldMask Mask, typename Func>
		AP_FORCE_INLINE  void for_each_pair (Func && f) const {
			// Iterate Chunks
			for (size_t c1 = range.start; c1 < range.end; ++c1) {

				// self-interaction (triangle loop within chunk)
				for (size_t i = 0; i < Container::chunk_size; ++i) {
					auto p1 = container.template restricted_at<Mask>(c1, i);
					for (size_t j = i + 1; j < Container::chunk_size; ++j) {
						auto p2 = container.template restricted_at<Mask>(c1, j);
						f(p1, p2);
					}
				}

				// pair-interaction (c1 with c2)
				for (size_t c2 = c1 + 1; c2 < range.end; ++c2) {
					for (size_t i = 0; i < Container::chunk_size; ++i) {
						auto p1 = container.template restricted_at<Mask>(c1, i);
						for (size_t j = 0; j < Container::chunk_size; ++j) {
							auto p2 = container.template restricted_at<Mask>(c2, j);
							f(p1, p2);
						}
					}
				}
			}
		}

		// Range represents chunk indices! (e.g., 0 to 4 means Chunks 0,1,2,3)
		Range range;
	private:
		Container & container;
	};
}



