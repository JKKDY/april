#pragma once

namespace april {
    enum class VectorPolicy {
    	Scalar,
    	Vector,	// force fully vectorized execution
    	Auto	// best effort vectorization
    };

    enum class ParallelPolicy {
    	Serial,
    	Threaded, // On-Node (Shared Memory) - e.g. std::thread, OpenMP, TBB
    	Distributed, // Multi-Node (Distributed Memory) - e.g. MPI, GASPI
    	Hybrid // Both Threaded and Distributed
    };

	template <VectorPolicy P>
	concept ExplicitVectorPolicy = (P == VectorPolicy::Scalar || P == VectorPolicy::Vector);

}




