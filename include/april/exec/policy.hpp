//  SPDX-License-Identifier: AGPL-3.0-only WITH ${PROJECT_NAME}-Commercial-Use-Exception-1.0
//  Copyright (c) ${YEAR} Julian Deller-Yee

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
		// what execution paths does a dispatcher support?
		enum class VectorTrait : uint8_t {
			ScalarOnly = 1 << 0,	// scalar path only
			VectorOnly = 1 << 1,	// vector path only
			Mixed = 1 << 2			// mixed execution path (e.g. peel + SIMD body)
		};

		AP_ENABLE_BITMASK_OPERATORS(VectorTrait)

		enum class ParallelTrait : uint8_t {
			None       = 0,			// no safe parallelism -> use atomic updates
			IntraBatch = 1 << 0,	// safe parallelism within a batch
			InterBatch = 1 << 1		// safe parallelism across batches
		};

		AP_ENABLE_BITMASK_OPERATORS(ParallelTrait)

		enum class ExecutionMode : uint8_t {
			Scalar = 1 << 0, // Capable of Scalar
			Vector = 1 << 1, // Capable of Vector
			Hybrid = Scalar | Vector // Generic lambda / Fully authorized
		};
		AP_ENABLE_BITMASK_OPERATORS(ExecutionMode)
	}

}




