#pragma once
#include "april/particle/particle_types.hpp"
#include "april/containers/layout/soa.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/batching/scalar_batch.hpp"
#include "april/containers/direct_sum/ds_batching.hpp"

namespace april::container::internal {

    template <class Config, class U>
    class DirectSumSoAImpl : public DirectSumCore<layout::SoA<Config, U>> {
    public:
        using Base = DirectSumCore<layout::SoA<Config, U>>;
        using SymmetricBatch = SymmetricParallelBatch<
            DirectSumSoAImpl,
            batching::AsymmetricScalarBatch<DirectSumSoAImpl>,
            batching::SymmetricScalarBatch<DirectSumSoAImpl>,
            exec::Executor
        >;
        using AsymmetricBatch =
            AsymmetricParallelBatch<
            DirectSumSoAImpl,
            batching::AsymmetricScalarBatch<DirectSumSoAImpl>,
            exec::Executor
        >;

        using Base::Base;
        friend Base;

        void generate_batches() {
            const auto n_types = static_cast<ParticleType>(this->bin_starts.size());

            // create batches for interacting particles of the same type
            for (ParticleType type = 0; type < n_types; type++) {
                auto [start, end, step] = this->get_physical_bin_range(type);
                if (end - start <= 1) continue;

                SymmetricBatch batch (*this, this->thread_executor);
                batch.types = {type, type};
                // batch.range = {start, end};
                batch.set_range({start,end});
                symmetric_batches.push_back(batch);
            }

            // create batches for interacting particles of different types
            for (ParticleType t1 = 0; t1 < n_types; t1++) {
                for (ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    auto [start1, end1, step1] = this->get_physical_bin_range(t1);
                    auto [start2, end2, step2] = this->get_physical_bin_range(t2);
                    if (start1 == end1 || start2 == end2) continue;

                    AsymmetricBatch batch (*this, this->thread_executor);
                    batch.types = {t1, t2};
                    // batch.range1 = {start1, end1};
                    // batch.range2 = {start2, end2};
                    batch.set_range({start1, end1}, {start2, end2});
                    asymmetric_batches.push_back(batch);
                }
            }
        }

        std::vector<SymmetricBatch> symmetric_batches;
        std::vector<AsymmetricBatch> asymmetric_batches;
    };
}

namespace april::container {
    struct DirectSumSoA {
        using ConfigT = DirectSumSoA;

        template <class U>
        using impl = internal::DirectSumSoAImpl<ConfigT, U>;
    };
}














