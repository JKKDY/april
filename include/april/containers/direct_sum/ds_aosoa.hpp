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

        using SymmetricBatch = batching::SymmetricChunkedBatch<DirectSumAoSoAImpl,  typename Base::ChunkT>;
        using AsymmetricBatch = batching::AsymmetricChunkedBatch<DirectSumAoSoAImpl, typename Base::ChunkT>;

        using SymTaskGroup = SymTaskGroup<SymmetricBatch, AsymmetricBatch>;
        using AsymTaskGroup = AsymTaskGroup<AsymmetricBatch>;


        using Base::Base;
        friend Base;

    private:
        using Base::bin_starts;
        using Base::bin_sizes;
        using Base::chunk_shift; // log_2(chunk_size) e.g. chunk_size = 8 -> chunk_shift = 0b100
        using Base::chunk_mask; // chunk_size -1 e.g. chunk_size = 8 -> chunk_mask = 0b111

        void generate_symmetric_group(this auto && self, ParticleType type) {
            using Derived = std::remove_cvref_t<decltype(self)>;
            using SymGroup = Derived::SymTaskGroup;

            auto range = self.get_physical_bin_range(type);
            if (range.size() <= 1) return;

            exec::BlockConfig config;
            config.alignment = self.chunk_size;

            auto schedule = exec::make_symmetric_schedule(range, config);
            SymGroup group;

            group.diagonals.reserve(schedule.diagonals.size());
            for (const auto& r : schedule.diagonals) {
                group.diagonals.push_back(self.create_symmetric_batch(type, r));
            }

            group.off_diagonals.resize(schedule.off_diagonals.size());
            for (size_t phase = 0; phase < schedule.off_diagonals.size(); ++phase) {
                for (const auto& [r1, r2] : schedule.off_diagonals[phase]) {
                    group.off_diagonals[phase].push_back(self.create_asymmetric_batch(type, r1, type, r2));
                }
            }

            self.sym_groups.push_back(std::move(group));
        }

        void generate_asymmetric_group(this auto && self, ParticleType type1, ParticleType type2) {
            using Derived = std::remove_cvref_t<decltype(self)>;
            using AsymGroup = Derived::AsymTaskGroup;

            auto range1 = self.get_physical_bin_range(type1);
            auto range2 = self.get_physical_bin_range(type2);
            if (range1.empty() || range2.empty()) return;

            exec::BlockConfig config;
            config.alignment = self.chunk_size;

            auto schedule = exec::make_bipartite_schedule(range1, range2, config);
            AsymGroup group;
            group.phases.resize(schedule.phases.size());

            for (size_t p = 0; p < schedule.phases.size(); ++p) {
                for (const auto& [r1, r2] : schedule.phases[p]) {
                    group.phases[p].push_back(self.create_asymmetric_batch(type1, r1, type2, r2));
                }
            }
            self.asym_groups.push_back(std::move(group));
        }

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
        using ConfigT = DirectSumAoSoA;

        template <class U>
        using impl = internal::DirectSumAoSoAImpl<ConfigT, U, ChunkSize>;
    };
}
