#pragma once
#include "april/particle/properties.hpp"
#include "april/containers/layout/aos.hpp"
#include "april/containers/direct_sum/ds_core.hpp"
#include "april/containers/batching/scalar_batch.hpp"
#include "april/containers/direct_sum/ds_batching.hpp"

namespace april::container::internal {

    template <class Config>
    class DirectSumAoSImpl : public DirectSumCore<layout::AoS<Config>> {
    public:
        using Base = DirectSumCore<layout::AoS<Config>>;

        using SymmetricBatch = batching::SymmetricScalarBatch<DirectSumAoSImpl, exec::ExecutionTrait::ScalarPath>;
        using AsymmetricBatch = batching::AsymmetricScalarBatch<DirectSumAoSImpl, exec::ExecutionTrait::ScalarPath>;

        using SymTaskGroup = SymmetricTaskGroup<SymmetricBatch, AsymmetricBatch>;
        using AsymTaskGroup = AsymmetricTaskGroup<AsymmetricBatch>;

        using Base::Base;
        friend Base;

    private:
        auto create_symmetric_batch(ParticleType type, const math::Range & range) {
            SymmetricBatch batch(*this);
            batch.types = {type, type};
            batch.range = range;
            return batch;
        }

        auto create_asymmetric_batch(ParticleType type1, const math::Range & range1, ParticleType type2, const math::Range & range2) {
            AsymmetricBatch batch(*this);
            batch.types = {type1, type2};
            batch.range1 = range1;
            batch.range2 = range2;
            return batch;
        }

        std::vector<SymTaskGroup> sym_groups;
        std::vector<AsymTaskGroup> asym_groups;
    };
}


namespace april::container {
    struct DirectSumAoS {
        template<class Config>
        using impl = internal::DirectSumAoSImpl<Config>;
    };
}














