#pragma once
#include "april/particle/particle_types.hpp"
#include "april/containers/layout/aosoa.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/direct_sum/ds_batching.hpp"
#include "april/containers/batching/chunked_batch.hpp"

namespace april::container::internal {
    template <class Config, size_t ChunkSize>
    class DirectSumAoSoAImpl : public DirectSumCore<layout::AoSoA<Config, ChunkSize>> {
    public:
        using Base = DirectSumCore<layout::AoSoA<Config, ChunkSize>>;
        using Base::chunk_size;

        using SymmetricBatch = batching::SymmetricChunkedBatch<DirectSumAoSoAImpl,  typename Base::ChunkT>;
        using AsymmetricBatch = batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>;

        using SymTaskGroup = SymmetricTaskGroup<SymmetricBatch, AsymmetricBatch>;
        using AsymTaskGroup = AsymmetricTaskGroup<AsymmetricBatch>;

        using Base::Base;
        friend Base;

    private:
        using Base::bin_starts;
        using Base::bin_sizes;
        using Base::chunk_shift; // log_2(chunk_size) e.g. chunk_size = 8 -> chunk_shift = 0b100
        using Base::chunk_mask; // chunk_size -1 e.g. chunk_size = 8 -> chunk_mask = 0b111

        auto create_symmetric_batch(ParticleType type, const math::Range& r) {
            const size_t c_start = r.start >> chunk_shift;
            const size_t c_stop  = (r.stop + chunk_size - 1) >> chunk_shift;

            SymmetricBatch batch(*this, this->ptr_chunks);
            batch.types = {type, type};
            batch.range_chunks = { c_start, c_stop };
            batch.range_tail = r.stop & chunk_mask;

            return batch;
        }

        auto create_asymmetric_batch(ParticleType t1, const math::Range& r1, ParticleType t2, const math::Range& r2) {
            AsymmetricBatch batch(*this, this->ptr_chunks);
            batch.types = {t1, t2};
            batch.range1_chunks = { r1.start >> chunk_shift, (r1.stop + chunk_size - 1) >> chunk_shift };
            batch.range2_chunks = { r2.start >> chunk_shift, (r2.stop + chunk_size - 1) >> chunk_shift };
            batch.range1_tail = r1.stop & chunk_mask;
            batch.range2_tail = r2.stop & chunk_mask;

            return batch;
        }

        std::vector<SymTaskGroup> sym_groups;
        std::vector<AsymTaskGroup> asym_groups;
    };
}


namespace april::container {
    template <size_t ChunkSize>
    struct DirectSumAoSoA {

        template <class Config>
        using impl = internal::DirectSumAoSoAImpl<Config, ChunkSize>;
    };
}
