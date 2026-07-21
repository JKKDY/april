#pragma once

#include <cstdint>
#include <utility>

#include "april/base/bitmask.hpp"
#include "april/utility/array.hpp"

namespace april {

    enum class VectorPolicy : std::uint8_t {
    	Scalar,
    	Vector,	// force fully vectorized execution
    	Auto	// best effort vectorization. Default.
    };

    enum class ParallelPolicy : std::uint8_t {
    	Serial,
    	Threaded
    };

	// future stub for distributed memory parallelism
	enum class DistributionPolicy : std::uint8_t {
		SingleProcess,
		Distributed
	};

	// future stub for accelerators e.g. GPUs
	enum class AcceleratorPolicy : std::uint8_t {
		Disabled,
		Preferred,
		Required
	};




	namespace exec {
		/**
		* @brief Kernel invocation modes.
		*
		* Values may be combined when one execution path invokes the kernel
		* using multiple modes.
		*/
		enum class ExecutionMode : uint8_t {
			None = 0,
			Scalar = 1 << 0,    // supports scalar execution
			Packed = 1 << 1,    // supports vector execution
			Group = 1 << 2		// Future stub: supports grouped execution (e.g. GPU SIMT)
		};
		APRIL_ENABLE_BITMASK_OPERATORS(ExecutionMode)

		/**
		* @brief Independently selectable execution paths offered by a dispatcher.
		*
		* Each template argument is the complete mode mask required by one path.
		* Paths are ordered by preference.
		*/
		template<ExecutionMode First, ExecutionMode... Rest>
		struct ExecutionPaths {
			static constexpr std::size_t size = 1 + sizeof...(Rest);
			static constexpr ExecutionMode first = First;
			static constexpr std::array<ExecutionMode, size> paths {First, Rest...};
		};

		namespace internal {
			template<typename T>struct is_execution_paths : std::false_type {};

			template<ExecutionMode... Paths>
			struct is_execution_paths<ExecutionPaths<Paths...>> : std::true_type {};


			// Create ExecutionPaths from an array
			template<class Left, class Right>
			consteval auto execution_paths_intersection() {
				constexpr auto common = utility::array_intersection<Left::paths, Right::paths>();
				static_assert(!common.empty(), "ExecutionPaths cannot represent an empty intersection");

				return [common]<std::size_t... I>(std::index_sequence<I...> ) {
					return ExecutionPaths<common[I]...>{};
				}(std::make_index_sequence<common.size()>{});
			}
		} // namespace internal

		template<typename T>
		concept IsExecutionPaths = internal::is_execution_paths<std::remove_cvref_t<T>>::value;

		template<class Left, class Right>
		using execution_paths_intersection_t = decltype(internal::execution_paths_intersection<Left, Right>());
	} // namespace exec






	namespace exec::internal {
	    /**
	     * @brief Checks whether all required modes are available.
	     */
	    [[nodiscard]]
	    constexpr bool supports_all(const ExecutionMode available, const ExecutionMode required) noexcept {
	        return (available & required) == required;
	    }

	    /**
	     * @brief Restricts kernel-supported modes according to VectorPolicy.
	     *
	     * VectorPolicy controls scalar versus packed host execution. Group mode is
	     * preserved because accelerator/SIMT selection is a separate policy axis.
	     */
	    template<VectorPolicy Policy, ExecutionMode SupportedModes >
	    consteval ExecutionMode allowed_execution_modes() {
	        if constexpr (Policy == VectorPolicy::Scalar) {
	            return SupportedModes & (ExecutionMode::Scalar | ExecutionMode::Group);
	        }
	        else if constexpr (Policy == VectorPolicy::Vector) {
	            return SupportedModes & (ExecutionMode::Packed | ExecutionMode::Group);
	        }
	        else {
	            return SupportedModes;
	        }
	    }

	    /**
	     * @brief Selects the first compatible execution path.
	     *
	     * A path is compatible when every mode it requires is present in
	     * AllowedModes. The order in ExecutionPaths defines path preference.
	     */
		template<ExecutionMode AllowedModes, IsExecutionPaths Paths>
		consteval ExecutionMode resolve_execution_path() {
	    	constexpr ExecutionMode selected = [] {
	    		for (const ExecutionMode path : Paths::paths) {
	    			if (supports_all(AllowedModes, path)) {
	    				return path;
	    			}
	    		}

	    		return ExecutionMode::None;
	    	}();

	    	static_assert(
				selected != ExecutionMode::None,
				"[APRIL] No compatible execution path exists between the "
				"batch, kernel, and selected execution policies."
			);

	    	return selected;
	    }
	} // namespace exec::internal
} // namespace april















