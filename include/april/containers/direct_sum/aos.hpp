#pragma once
#include "april/particle/defs.hpp"
#include "april/containers/layout/aos.hpp"
#include "april/containers/direct_sum/core.hpp"
#include "april/containers/batching.hpp"

namespace april::container::direct_sum {

    template <class Config, class U>
    class DirectSumAoSImpl : public DirectSumCore<layout::AoS<Config, U>> {
    public:
        using Base = DirectSumCore<layout::AoS<Config, U>>;
        using SymmetricBatch = batching::SymmetricScalarBatch<DirectSumAoSImpl>;
        using AsymmetricBatch = batching::AsymmetricScalarBatch<DirectSumAoSImpl>;

        using Base::Base;
        friend Base;

        void generate_batches(const std::vector<std::pair<size_t, size_t>> & type_ranges) {
            const auto n_types = static_cast<env::ParticleType>(type_ranges.size());
            for (env::ParticleType type = 0; type < n_types; type++) {
                SymmetricBatch batch (*this);
                batch.types = {type, type};
                batch.range = type_ranges[type];
                symmetric_batches.push_back(batch);
            }

            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    AsymmetricBatch batch (*this);
                    batch.types = {t1, t2};
                    batch.range1 = type_ranges[t1];
                    batch.range2 = type_ranges[t2];
                    asymmetric_batches.push_back(batch);
                }
            }
        }

        std::vector<SymmetricBatch> symmetric_batches;
        std::vector<AsymmetricBatch> asymmetric_batches;
    };

    struct DirectSumAoS {
        using ConfigT = DirectSumAoS;

        template <class U>
        using impl = DirectSumAoSImpl<ConfigT, U>;
    };


}