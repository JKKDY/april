#pragma once
#include <concepts>

#include "april/exec/policy.hpp"
#include "april/particle/properties.hpp"


namespace april::container::batching {


	//-----------------------
	// CONVENIENCE DEFINITION
	//-----------------------
	template<
		std::size_t Arity, // future stub for tri-wise or higher order interactions. Currently only arity=2 is in use
		exec::IsExecutionPaths Paths
	>
	struct BatchBase {
		static_assert(Arity > 0, "A batch must have a nonzero arity.");

		static constexpr std::size_t arity = Arity;

		using execution_paths = Paths;

		std::array<ParticleType, Arity> types{};
	};


	//--------------
	// BATCH CONCEPT
	//--------------
	namespace internal {
		template<typename B>
		concept IsBatchImpl = requires(const B& batch) {
			// Must declare a valid execution-path set.
			typename B::execution_paths;
			requires exec::IsExecutionPaths<typename B::execution_paths>;

			// Must declare its arity.
			{ B::arity } -> std::convertible_to<std::size_t>;

			// Must expose one ParticleType per argument.
			{ batch.types } -> std::convertible_to<const std::array<ParticleType, B::arity>&>;

			// Must implement at least the first declared execution path.
			batch.template for_each<B::execution_paths::first>(
				universal_kernel([](auto&&...) {})
			);
		};

		template<typename B>
		concept IsTopologyBatchImpl = requires(const B& batch) {
			typename B::execution_paths;
			requires exec::IsExecutionPaths<typename B::execution_paths>;

			// Must declare its arity.
			{ B::arity } -> std::convertible_to<std::size_t>;

			// Must expose one ParticleType per argument.
			{ batch.representatives } -> std::convertible_to<const std::array<ParticleType, B::arity>&>;
			{ batch.interactions } -> std::convertible_to<const std::vector<std::array<ParticleID, B::arity>>&>;

			// Must implement at least the first declared execution path.
			batch.template for_each<B::execution_paths::first>(
				universal_kernel([](auto&&...) {})
			);
		};
	} // namespace internal

	template<typename T>
	concept IsBatch = internal::IsBatchImpl<std::remove_cvref_t<T>>;;

	template<typename T>
	concept IsTopologyBatch = internal::IsTopologyBatchImpl<std::remove_cvref_t<T>>;


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






