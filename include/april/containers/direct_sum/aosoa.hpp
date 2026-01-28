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
        using SymmetricBatch = batching::SymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::Chunk>;
        using AsymmetricBatch = batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::Chunk>;

        using Base::Base;
        friend Base;

        void generate_batches(const std::vector<std::pair<size_t, size_t>> &) {
            const auto n_types = static_cast<env::ParticleType>(this->bin_starts.size() - 1); // bin_starts has N+1 entries

            for (env::ParticleType type = 0; type < n_types; type++) {
                auto [chunk_start, chunk_end] = this->get_chunk_bin_range(type);

                if (chunk_start == chunk_end) continue;

                SymmetricBatch batch (*this,  this->ptr_chunks);
                batch.types = {type, type};
                batch.range_chunks = {chunk_start, chunk_end};
                batch.range_tail = this->bin_sizes[type] % chunk_size;

                symmetric_batches.push_back(batch);
            }

            // build asymmetric batches
            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    auto [c1_start, c1_end] = this->get_chunk_bin_range(t1);
                    auto [c2_start, c2_end] = this->get_chunk_bin_range(t2);

                    if (c1_start == c1_end || c2_start == c2_end) continue;

                    AsymmetricBatch batch (*this, this->ptr_chunks);
                    batch.types = {t1, t2};
                    batch.range1_chunks = {c1_start, c1_end};
                    batch.range2_chunks = {c2_start, c2_end};
                    batch.range1_tail = this->bin_sizes[t1] % chunk_size;
                    batch.range2_tail = this->bin_sizes[t2] % chunk_size;

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