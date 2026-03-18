#pragma once

#include<vector>

#include "april/exec/policy.hpp"
#include "april/particle/particle_types.hpp"


namespace april::container::batching {
    template <typename Container>
    struct TopologyBatch {
        static constexpr auto parallel_trait = exec::ParallelTrait::None;
        static constexpr auto vector_trait = exec::VectorTrait::ScalarPath;

        std::pair<ParticleType, ParticleType> representatives;
        std::vector<std::pair<ParticleID, ParticleID>> pairs;
        Container* container_ptr;

        template<ParallelPolicy P = ParallelPolicy::Serial, typename Kernel>
        void for_each_pair(Kernel&& kernel) const {
            for (const auto& [id1, id2] : pairs) {
                auto p1 = container_ptr->template at_id<Kernel::Read, Kernel::Write>(id1);
                auto p2 = container_ptr->template at_id<Kernel::Read, Kernel::Write>(id2);
                kernel(p1, p2);
            }
        }
    };
}
