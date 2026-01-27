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
        using SymmetricBatch = batching::SymmetricChunkedBatch<DirectSumAoSoAImpl>;
        using AsymmetricBatch = batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl>;

        using Base::Base;
        friend Base;

        void generate_batches(const std::vector<std::pair<size_t, size_t>> & /* unused_particle_ranges */) {
            // IGNORE the input ranges! They are particle-based.
            // We need the CHUNK ranges that we just built in reorder_storage.

            // Assume we can access 'bin_starts' or a helper from the base class AoSoA.
            // Since DirectSumAoSoAImpl inherits DirectSumCore -> AoSoA, we can call 'this->bin_range(type)'.

            const auto n_types = static_cast<env::ParticleType>(this->bin_starts.size() - 1); // bin_starts has N+1 entries

            for (env::ParticleType type = 0; type < n_types; type++) {
                auto [chunk_start, chunk_end] = this->bin_range(type);

                if (chunk_start == chunk_end) continue; // Empty type

                SymmetricBatch batch (*this);
                batch.types = {type, type};
                batch.range = {chunk_start, chunk_end};
                symmetric_batches.push_back(batch);
            }

            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    auto [c1_start, c1_end] = this->bin_range(t1);
                    auto [c2_start, c2_end] = this->bin_range(t2);

                    if (c1_start == c1_end || c2_start == c2_end) continue;

                    AsymmetricBatch batch (*this);
                    batch.types = {t1, t2};
                    batch.range1 = {c1_start, c1_end};
                    batch.range2 = {c2_start, c2_end};
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