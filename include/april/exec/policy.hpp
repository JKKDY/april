#pragma once

#include <cstdint>
#include <utility>

#include "april/base/bitmask.hpp"

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
		// describes the structural execution path(s) offered by the dispatcher
		enum class ExecutionTrait : uint8_t {
			ScalarPath = 1 << 0,  // dispatcher provides a scalar execution path -> kernel must support Scalar
			VectorPath = 1 << 1,  // dispatcher provides a vector execution path -> kernel must support Vector
			HybridPath = 1 << 2   // scalar + vector occur in the same path (e.g. SIMD body + scalar tail) -> kernel must support both
		};
		APRIL_ENABLE_BITMASK_OPERATORS(ExecutionTrait)

		// describes which execution modes a kernel implementation supports
		enum class ExecutionMode : uint8_t {
			Scalar = 1 << 0,    // supports scalar execution
			Packed = 1 << 1,    // supports vector execution
			Group = 1 << 2		// Future stub: supports grouped execution (e.g. GPU SIMT)
		};
		APRIL_ENABLE_BITMASK_OPERATORS(ExecutionMode)


	}


	namespace exec::internal {
		// given a policy and execution mode(s), return a valid mode
		template <VectorPolicy Policy, ExecutionMode Mode>
		constexpr ExecutionMode allowed_execution_modes() {
			// user requests vectorized execution
			if constexpr (Policy == VectorPolicy::Vector) {
				static_assert(static_cast<bool>(Mode & ExecutionMode::Packed),
					"Error: VectorPolicy::Vector was requested, but the provided Kernel does not support Vector execution.");
				return ExecutionMode::Packed;
			}
			// user requests scalar execution
			else if constexpr (Policy == VectorPolicy::Scalar) {
				static_assert(static_cast<bool>(Mode & ExecutionMode::Scalar),
					"Error: VectorPolicy::Scalar was requested, but the provided Kernel does not support Scalar execution.");
				return ExecutionMode::Scalar;
			}
			// auto mode
			else {
				return Mode;
			}
		}

		template <ExecutionTrait Trait>
		consteval ExecutionMode required_execution_modes() {
			if constexpr (static_cast<bool>(Trait & ExecutionTrait::HybridPath)) {
				return ExecutionMode::Scalar | ExecutionMode::Packed; // requires both scalar and vector
			}
			if constexpr (static_cast<bool>(Trait & ExecutionTrait::ScalarPath) &&
				static_cast<bool>(Trait & ExecutionTrait::VectorPath)) {
				return  static_cast<ExecutionMode>(0); // dont care which mode is used
			}
			if constexpr (static_cast<bool>(Trait & ExecutionTrait::ScalarPath)) {
				return ExecutionMode::Scalar;
			}
			if constexpr (static_cast<bool>(Trait & ExecutionTrait::VectorPath)) {
				return ExecutionMode::Packed;
			}
			std::unreachable();
		}

		template <ExecutionMode Valid, ExecutionMode Required>
		consteval ExecutionMode resolve_execution_mode() {
			// Batch has strict structural demands. We must follow them.
			if constexpr (Required != static_cast<ExecutionMode>(0)) {
				return Required;
			}

			// Batch provides independent paths (Required == 0).
			if constexpr (static_cast<bool>(Valid & ExecutionMode::Packed)) {
				return ExecutionMode::Packed; // Prefer Vector
			} else {
				return ExecutionMode::Scalar; // Fallback to pure Scalar
			}
		}
	}

}















