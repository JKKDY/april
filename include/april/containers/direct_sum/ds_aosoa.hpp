#pragma once
#include "april/particle/defs.hpp"
#include "april/containers/layout/aosoa.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/batching/chunked.hpp"

namespace april::container::internal {

    template <class Config, class U>
    class DirectSumAoSoAImpl : public DirectSumCore<layout::AoSoA<Config, U>> {
    public:
        using Base = DirectSumCore<layout::AoSoA<Config, U>>;
        using Base::chunk_size;
        using SymmetricBatch = SymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::Chunk>;
        using AsymmetricBatch = AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::Chunk>;

        using Base::Base;
        friend Base;

        void generate_batches() {
            const auto n_types = static_cast<env::ParticleType>(this->bin_starts.size() - 1); // bin_starts has N+1 entries

            for (env::ParticleType type = 0; type < n_types; type++) {
                auto range_chunks = this->get_chunk_bin_range(type);

                if (range_chunks.empty()) continue;

                SymmetricBatch batch (*this,  this->ptr_chunks);
                batch.types = {type, type};
                batch.range_chunks = range_chunks;
                batch.range_tail = this->bin_sizes[type] % chunk_size;

                symmetric_batches.push_back(batch);
            }

            // build asymmetric batches
            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    auto range_chunks1 = this->get_chunk_bin_range(t1);
                    auto range_chunks2 = this->get_chunk_bin_range(t2);

                    if (range_chunks1.empty() || range_chunks2.empty()) continue;

                    AsymmetricBatch batch (*this, this->ptr_chunks);
                    batch.types = {t1, t2};
                    batch.range1_chunks = range_chunks1;
                    batch.range2_chunks = range_chunks2;
                    batch.range1_tail = this->bin_sizes[t1] % chunk_size;
                    batch.range2_tail = this->bin_sizes[t2] % chunk_size;

                    asymmetric_batches.push_back(batch);
                }
            }
        }

        std::vector<SymmetricBatch> symmetric_batches;
        std::vector<AsymmetricBatch> asymmetric_batches;
    };
}


namespace april::container {
    struct DirectSumAoSoA {
        using ConfigT = DirectSumAoSoA;

        template <class U>
        using impl = internal::DirectSumAoSoAImpl<ConfigT, U>;
    };
}