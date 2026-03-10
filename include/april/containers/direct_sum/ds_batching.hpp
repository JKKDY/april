#pragma once

#include <vector>

#include "april/exec/policy.hpp"
#include "april/exec/particle_kernel.hpp"
#include "april/containers/batching/common.hpp"
#include "april/math/range.hpp"
#include "april/exec/concurrency.hpp"

namespace april::container::internal {

    template<typename Container, typename AsymmetricBatch, typename SymmetricBatch>
    struct SymmetricParalleldBatch : batching::BatchBase<exec::ParallelTrait::IntraBatch,
        AsymmetricBatch::vector_trait & SymmetricBatch::vector_trait>
    {
        explicit SymmetricParalleldBatch(Container & container) : container(container) {}

        template<ParallelPolicy P, exec::ExecutionMode E, exec::IsKernel Func>
        void for_each_pair (Func && f) const {
            for (auto & block : diagonal_phase) {
                block.template for_each_pair<P, E>(f);
            }

            for (auto & phase : off_diagonal_phases) {
                for (auto & block : phase) {
                    block.template for_each_pair<P, E>(f);
                }
            }
        }

        void set_range(const math::Range & range, const size_t oversubscription = 2) {
            const size_t n_threads = exec::CPU_THREADS;
            size_t B = std::max<size_t>(1, n_threads * oversubscription);
            if (B % 2 != 0) B++;

            // Partition the range into B blocks
            std::vector<math::Range> blocks(B);
            const size_t chunk_size = range.size() / B;
            const size_t remainder = range.size() % B;
            size_t current_idx = range.start;

            for (size_t i = 0; i < B; ++i) {
                const size_t size = chunk_size + (i < remainder ? 1 : 0);
                blocks[i] = {current_idx, current_idx + size};
                current_idx += size;
            }

            // Phase 0: create symmetric batches for the diagonals
            for (size_t i = 0; i < B; ++i) {
                if (blocks[i].start == blocks[i].stop) continue;

                SymmetricBatch sym(this->container);
                sym.types = this->types;
                sym.range = blocks[i];
                diagonal_phase.push_back(sym);
            }

            // Phases 1 to B-1: create asymmetric batches for off diagonals
            off_diagonal_phases.resize(B - 1);
            std::vector<size_t> circle(B);
            for (size_t i = 0; i < circle.size(); ++i) circle[i] = i;

            for (size_t phase = 0; phase < B - 1; ++phase) {
                for (size_t i = 0; i < B / 2; ++i) {
                    // block b1 interacts with block b2
                    size_t b1 = circle[i];
                    size_t b2 = circle[B - 1 - i];

                    const size_t idx1 = std::min(b1, b2);
                    const size_t idx2 = std::max(b1, b2);

                    // if empty skip
                    if (blocks[idx1].start == blocks[idx1].stop ||
                        blocks[idx2].start == blocks[idx2].stop) continue;

                    AsymmetricBatch asym(container);
                    asym.types = this->types;
                    asym.range1 = blocks[idx1];
                    asym.range2 = blocks[idx2];
                    off_diagonal_phases[phase].push_back(asym);
                }

                // rotate circle array for next phase: keep index 0 fixed, rotate 1..B-1 right
                const size_t last = circle.back();
                for (size_t i = B - 1; i > 1; --i) {
                    circle[i] = circle[i - 1];
                }
                circle[1] = last;
            }
        }

    private:
        Container & container;
        std::vector<SymmetricBatch> diagonal_phase;
        std::vector<std::vector<AsymmetricBatch>> off_diagonal_phases;
    };
}















