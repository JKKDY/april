#pragma once
#include "april/particle/defs.hpp"
#include "april/containers/layout/aosoa.hpp"
#include "april/containers/direct_sum/core.hpp"
#include "april/containers/batching.hpp"

namespace april::container::direct_sum {

    template <class Config, class U>
    class DirectSumAoSoAImpl : public DirectSumCore<layout::AoSoA<Config, U>> {
    public:
        using Base = DirectSumCore<layout::AoSoA<Config, U>>;
        using Base::chunk_size;
        using SymmetricBatch = batching::SymmetricScalarBatch<DirectSumAoSoAImpl>;
        using AsymmetricBatch = batching::AsymmetricScalarBatch<DirectSumAoSoAImpl>;

        using Base::Base;
        friend Base;

        void generate_batches(const std::vector<std::pair<size_t, size_t>> &) {
            const auto n_types = static_cast<env::ParticleType>(this->bin_counts.size());

            // build symmetric batches
            for (env::ParticleType type = 0; type < n_types; type++) {
                // Use the new helper to get PHYSICAL indices
                auto [start, end] = this->get_physical_bin_range(type);

                if (start == end) continue;

                SymmetricBatch batch (*this);
                batch.types = {type, type};
                batch.range = {start, end};
                symmetric_batches.push_back(batch);
            }

            // build symmetric batches
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

    struct DirectSumAoSoA {
        using ConfigT = DirectSumAoSoA;

        template <class U>
        using impl = DirectSumAoSoAImpl<ConfigT, U>;
    };


}