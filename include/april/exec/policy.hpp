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
		// what execution paths does a dispatcher have?
		enum class VectorTrait : uint8_t {
			ScalarOnly = 1 << 0,	// has scalar path only
			VectorOnly = 1 << 1,	// has vector path only
			Mixed = 1 << 2			// has mixed execution path (e.g. peel + SIMD body)
		};
		AP_ENABLE_BITMASK_OPERATORS(VectorTrait)

		// what parallelization paths does a dispatcher have?
		enum class ParallelTrait : uint8_t {
			None       = 0,			// no safe parallelism -> use atomic updates
			IntraBatch = 1 << 0,	// safe parallelism within a batch
			InterBatch = 1 << 1		// safe parallelism across batches
		};
		AP_ENABLE_BITMASK_OPERATORS(ParallelTrait)

		// what mode should can a kernel run as?
		enum class ExecutionMode : uint8_t {
			Scalar = 1 << 0, // Capable of Scalar
			Vector = 1 << 1, // Capable of Vector
			Hybrid = Scalar | Vector // Generic - can do both scalar and vector
		};
		AP_ENABLE_BITMASK_OPERATORS(ExecutionMode)


		template <VectorPolicy Policy, ExecutionMode Mode>
		constexpr ExecutionMode resolve_execution_mode() {
			if constexpr (Policy == VectorPolicy::Vector) {
				static_assert(static_cast<bool>(Mode & ExecutionMode::Vector),
					"Error: VectorPolicy::Vector was requested, but the provided Kernel does not support Vector execution.");
				return ExecutionMode::Vector;
			}
			else if constexpr (Policy == VectorPolicy::Scalar) {
				static_assert(static_cast<bool>(Mode & ExecutionMode::Scalar),
					"Error: VectorPolicy::Scalar was requested, but the provided Kernel does not support Scalar execution.");
				return ExecutionMode::Scalar;
			}
			else {
				return Mode;
			}
		}
	}

}













