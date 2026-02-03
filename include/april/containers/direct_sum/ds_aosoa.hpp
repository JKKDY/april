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
        using SymmetricBatch = SymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>;
        using AsymmetricBatch = AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>;

        using Base::Base;
        friend Base;

    private:
        using Base::bin_starts;
        using Base::bin_sizes;
        using Base::chunk_shift;
        using Base::chunk_mask;

        void generate_batches() {
            const auto n_types = static_cast<env::ParticleType>(bin_starts.size() - 1); // bin_starts has N+1 entries

            for (env::ParticleType type = 0; type < n_types; type++) {
                const size_t start = bin_starts[type];
                const size_t end = bin_starts[type+1];
                const size_t size = bin_sizes[type];

                if (size <= 1) continue;

                SymmetricBatch batch (*this,  this->ptr_chunks);
                batch.types = {type, type};
                batch.range_chunks = {start >> chunk_shift, end >> chunk_shift};
                batch.range_tail = size & chunk_mask;

                symmetric_batches.push_back(batch);
            }

            // build asymmetric batches
            for (env::ParticleType t1 = 0; t1 < n_types; t1++) {
                for (env::ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    const size_t start1 = bin_starts[t1];
                    const size_t size1 = bin_sizes[t1];
                    const size_t end1 = bin_starts[t1+1];

                    const size_t start2 = bin_starts[t2];
                    const size_t size2 = bin_sizes[t2];
                    const size_t end2 = bin_starts[t2+1];

                    if (size1==0 || size2==0) continue;

                    AsymmetricBatch batch (*this, this->ptr_chunks);
                    batch.types = {t1, t2};
                    batch.range1_chunks = {start1 >> chunk_shift, end1 >> chunk_shift};
                    batch.range2_chunks = {start2 >> chunk_shift, end2 >> chunk_shift};
                    batch.range1_tail = size1 & chunk_mask;
                    batch.range2_tail = size2 & chunk_mask;

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