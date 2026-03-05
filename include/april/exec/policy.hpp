#pragma once

#include <cstdint>
#include "april/base/bitmask.hpp"

namespace april {
    enum class VectorPolicy {
    	Scalar,
    	Vector,	// force fully vectorized execution
    	Auto	// best effort vectorization
    };

    enum class ParallelPolicy {
    	Serial,
    	Threaded, // On-Node (Shared Memory) - e.g. std::thread, OpenMP, TBB
    	Distributed, // Multi-Node (Distributed Memory) - e.g. MPI
    	Hybrid // Both Threaded and Distributed
    };


	namespace exec::internal {
		// describes the structural execution path(s) offered by the dispatcher
		enum class VectorTrait : uint8_t {
			ScalarPath = 1 << 0,   // dispatcher provides a scalar execution path -> kernel must support Scalar
			VectorPath = 1 << 1,   // dispatcher provides a vector execution path -> kernel must support Vector
			HybridPath = 1 << 2    // scalar + vector occur in the same path (e.g. peel + SIMD body) -> kernel must support both
		};
		AP_ENABLE_BITMASK_OPERATORS(VectorTrait)

		// describes which parallel execution strategies a dispatcher supports
		enum class ParallelTrait : uint8_t {
			None       = 0,			// no safe parallelism -> requires atomic updates
			IntraBatch = 1 << 0,	// safe parallelism within a batch
			InterBatch = 1 << 1		// safe parallelism across batches
		};
		AP_ENABLE_BITMASK_OPERATORS(ParallelTrait)

		// describes which execution modes a kernel implementation supports
		enum class ExecutionMode : uint8_t {
			Scalar = 1 << 0,       // supports scalar execution
			Vector = 1 << 1,       // supports vector execution
			Hybrid = Scalar | Vector // supports both scalar and vector execution
		};
		AP_ENABLE_BITMASK_OPERATORS(ExecutionMode)


		// given a policy and execution mode(s), return a valid mode
		template <VectorPolicy Policy, ExecutionMode Mode>
		constexpr ExecutionMode valid_execution_modes() {
			// user requests vectorized execution
			if constexpr (Policy == VectorPolicy::Vector) {
				static_assert(static_cast<bool>(Mode & ExecutionMode::Vector),
					"Error: VectorPolicy::Vector was requested, but the provided Kernel does not support Vector execution.");
				return ExecutionMode::Vector;
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

		template <VectorTrait Trait>
		consteval ExecutionMode required_execution_modes() {
			if constexpr (static_cast<bool>(Trait & VectorTrait::HybridPath)) {
				return ExecutionMode::Hybrid; // requires both scalar and vector
			}
			if constexpr (static_cast<bool>(Trait & VectorTrait::ScalarPath) &&
				static_cast<bool>(Trait & VectorTrait::VectorPath)) {
				return  static_cast<ExecutionMode>(0); // dont care which mode is used
			}
			if constexpr (static_cast<bool>(Trait & VectorTrait::ScalarPath)) {
				return ExecutionMode::Scalar;
			}
			if constexpr (static_cast<bool>(Trait & VectorTrait::VectorPath)) {
				return ExecutionMode::Vector;
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
			if constexpr (static_cast<bool>(Valid & ExecutionMode::Vector)) {
				return ExecutionMode::Vector; // Prefer Vector
			} else {
				return ExecutionMode::Scalar; // Fallback to pure Scalar
			}
		}
	}

}













