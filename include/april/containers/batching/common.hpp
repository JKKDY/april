#pragma once
#include <concepts>

#include "april/exec/policy.hpp"
#include "april/particle/properties.hpp"


namespace april::container::batching {


	//-----------------------
	// CONVENIENCE DEFINITION
	//-----------------------
	template<exec::ExecutionTrait V>
	struct BatchBase {
		static constexpr auto vector_trait = V;
		std::pair<ParticleType, ParticleType> types {};
	};


	//--------------
	// BATCH CONCEPT
	//--------------
	template<typename T>
	concept IsBatch =requires(const T& b) {
		// must have static constexpr trait flags
		{ T::vector_trait }	-> std::convertible_to<exec::ExecutionTrait>;

		// must have type pair
		{ b.types } -> std::convertible_to<std::pair<ParticleType, ParticleType>>;

		// must have vector trait exists
		{ std::remove_cvref_t<T>::vector_trait } -> std::convertible_to<exec::ExecutionTrait>;

		// must have a for_each_pair function
		b.template for_each_pair<exec::ExecutionMode::Scalar | exec::ExecutionMode::Packed>(
			universal_kernel([](auto&&, auto&&) {})
		);
	};

	template <typename T>
	concept IsTopologyBatch =
	requires(T t) {
		{ t.representatives } -> std::convertible_to<std::pair<ParticleType, ParticleType>>;

		t.template for_each_pair<ParallelPolicy::Serial>(
			scalar_kernel([](auto&&, auto&&) {})
		);
	};


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
			return v;
		}
	};
}






