#pragma once
#include "april/particle/particle_types.hpp"
#include "april/containers/layout/aosoa.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/batching/chunked_batch.hpp"

namespace april::container::internal {
    template <class Config, class U, size_t ChunkSize>
    class DirectSumAoSoAImpl : public DirectSumCore<layout::AoSoA<Config, U, ChunkSize>> {
    public:
        using Base = DirectSumCore<layout::AoSoA<Config, U, ChunkSize>>;
        using Base::chunk_size;

        using SymmetricBatch = SymmetricParallelChunkedBatch<
            DirectSumAoSoAImpl,
            batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>,
            batching::SymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>,
            typename Base::ChunkT,
            exec::Executor
        >;
        using AsymmetricBatch = AsymmetricParallelChunkedBatch<
            DirectSumAoSoAImpl,
            batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>,
            typename Base::ChunkT,
            exec::Executor
        >;


        using Base::Base;
        friend Base;

    private:
        using Base::bin_starts;
        using Base::bin_sizes;
        using Base::chunk_shift;
        using Base::chunk_mask;

        void generate_batches() {
            const auto n_types = static_cast<ParticleType>(bin_starts.size() - 1); // bin_starts has N+1 entries

            for (ParticleType type = 0; type < n_types; type++) {
                const size_t start = bin_starts[type];
                const size_t end = bin_starts[type + 1];
                const size_t size = bin_sizes[type];

                if (size <= 1) continue;

                math::Range range = {start >> chunk_shift, end >> chunk_shift};
                size_t tail = size & chunk_mask;

                SymmetricBatch batch(*this, this->ptr_chunks, this->thread_executor);
                batch.types = {type, type};

                batch.set_range(range, tail);

                symmetric_batches.push_back(batch);
            }

            // build asymmetric batches
            for (ParticleType t1 = 0; t1 < n_types; t1++) {
                for (ParticleType t2 = t1 + 1; t2 < n_types; t2++) {
                    const size_t start1 = bin_starts[t1];
                    const size_t size1 = bin_sizes[t1];
                    const size_t end1 = bin_starts[t1 + 1];

                    const size_t start2 = bin_starts[t2];
                    const size_t size2 = bin_sizes[t2];
                    const size_t end2 = bin_starts[t2 + 1];

                    if (size1 == 0 || size2 == 0) continue;

                    const math::Range range1 = {start1 >> chunk_shift, end1 >> chunk_shift};
                    const size_t tail1 = size1 & chunk_mask;

                    const math::Range range2 = {start2 >> chunk_shift, end2 >> chunk_shift};
                    const size_t tail2 = size2 & chunk_mask;

                    AsymmetricBatch batch(*this, this->ptr_chunks, this->thread_executor);
                    batch.types = {t1, t2};
                    batch.set_range(range1, tail1, range2, tail2);

                    asymmetric_batches.push_back(batch);
                }
            }
        }

        std::vector<SymmetricBatch> symmetric_batches;
        std::vector<AsymmetricBatch> asymmetric_batches;
    };
}


namespace april::container {
    template <size_t ChunkSize>
    struct DirectSumAoSoA {
        using ConfigT = DirectSumAoSoA;

        template <class U>
        using impl = internal::DirectSumAoSoAImpl<ConfigT, U, ChunkSize>;
    };
}
