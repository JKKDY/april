#pragma once
#include "april/particle/defs.hpp"
#include "april/containers/layout/aos.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/batching/scalar.hpp"

namespace april::container::internal {

    template <class Config, class U>
    class DirectSumAoSImpl : public DirectSumCore<layout::AoS<Config, U>> {
    public:
        using Base = DirectSumCore<layout::AoS<Config, U>>;
        using SymmetricBatch = SymmetricScalarBatch<DirectSumAoSImpl>;
        using AsymmetricBatch = AsymmetricScalarBatch<DirectSumAoSImpl>;

        using Base::Base;
        friend Base;

        void generate_batches() {
            const auto n_types = static_cast<env::ParticleType>(this->bin_starts.size());
            for (env::ParticleType type = 0; type < n_types; type++) {
                auto [start, end] = this->get_physical_bin_range(type);
                if (end - start <= 1) continue;

                SymmetricBatch batch (*this);
                batch.types = {type, type};
                batch.range = {start, end};
                symmetric_batches.push_back(batch);
            }

            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    auto [start1, end1] = this->get_physical_bin_range(t1);
                    auto [start2, end2] = this->get_physical_bin_range(t2);
                    if (start1 == end1 || start2 == end2) continue;

                    AsymmetricBatch batch (*this);
                    batch.types = {t1, t2};
                    batch.range1 = {start1, end1};
                    batch.range2 = {start2, end2};
                    asymmetric_batches.push_back(batch);
                }
            }
        }

        std::vector<SymmetricBatch> symmetric_batches;
        std::vector<AsymmetricBatch> asymmetric_batches;
    };
}


namespace april::container {
    struct DirectSumAoS {
        using ConfigT = DirectSumAoS;

        template <class U>
        using impl = internal::DirectSumAoSImpl<ConfigT, U>;
    };
}