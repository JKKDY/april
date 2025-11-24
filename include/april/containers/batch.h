#pragma once
#include <cstdint>
#include <utility>

#include "april/common.h"
#include "april/particle/defs.h"


namespace april::container {
	enum class BatchSymmetry : uint8_t {
		symmetric,
		asymmetric
	};

	enum class BatchRegion : uint8_t {
		boundary,
		inner
	};

	// Base constraints common to all batches
	template <typename T>
	concept IsBatchBase = requires(const T& b) {
		// Must have static constexpr configuration flags
		{ T::symmetry } -> std::convertible_to<BatchSymmetry>;
		{ T::region } -> std::convertible_to<BatchRegion>;
		{ T::parallelize_inner } -> std::convertible_to<bool>;

		// Must have type pair
		{ b.types } -> std::convertible_to<std::pair<env::ParticleType, env::ParticleType>>;
	};

	// Concept for Symmetric Batches
	// Enforces single index span and symmetry flag
	template <typename T>
	concept IsSymmetricBatch = IsBatchBase<T> &&
		(T::symmetry == BatchSymmetry::symmetric) &&
		requires (const T& b) {
		std::ranges::input_range<decltype(b.type_indices)>; // must be iterable
		std::same_as<std::ranges::range_value_t<decltype(b.type_indices)>, size_t>; // iteration yields size_t
		};

	// Concept for Asymmetric Batches
	// Enforces two distinct index spans and symmetry flag
	template <typename T>
	concept IsAsymmetricBatch = IsBatchBase<T> &&
		(T::symmetry == BatchSymmetry::asymmetric) &&
		requires (const T& b) {
		std::ranges::input_range<decltype(b.type1_indices)>;
		std::same_as<std::ranges::range_value_t<decltype(b.type1_indices)>, size_t>;

		std::ranges::input_range<decltype(b.type2_indices)>;
		std::same_as<std::ranges::range_value_t<decltype(b.type2_indices)>, size_t>;
		};

	template <typename T>
	concept IsBatch = IsSymmetricBatch<T> || IsAsymmetricBatch<T>;

	struct NoBatchBCP {
		template <class T>
		constexpr T operator()(T&& v) const noexcept {
			return std::forward<T>(v); // identity
		}
	};

	template<typename F>
	concept IsBCP = requires(const F& f, const vec3& v) {
		{ f(v) } -> std::convertible_to<vec3>;
	};
}
