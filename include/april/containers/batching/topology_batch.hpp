#pragma once

#include <vector>
#include <utility>
#include <numeric>

#include "april/exec/policy.hpp"
#include "april/particle/properties.hpp"
#include "april/utility/graph.hpp"
#include "april/exec/kernel.hpp"

namespace april::container::batching {

    using ParticleGraph = utility::graph::EdgeList<ParticleID>;
    using ParticleSet = std::vector<ParticleID>;


    template<std::size_t Arity, typename Container>
    struct TopologyBatch {
        static constexpr std::size_t arity = Arity;

        using execution_paths = exec::ExecutionPaths<exec::ExecutionMode::Scalar>;

        std::array<ParticleType, Arity> representatives{};
        std::vector<std::array<ParticleID, Arity>> interactions;
        Container* container_ptr = nullptr;

        template<exec::ExecutionMode Mode, exec::IsKernel Kernel>
        void for_each(Kernel&& kernel) const {
            static_assert(Mode == exec::ExecutionMode::Scalar, "TopologyBatch only supports scalar execution.");

            for (const auto& ids : interactions) {
                invoke(
                    ids,
                    std::forward<Kernel>(kernel),
                    std::make_index_sequence<Arity>{}
                );
            }
        }

    private:
        template<exec::IsKernel Kernel, std::size_t... I>
        void invoke(
            const std::array<ParticleID, Arity>& ids,
            Kernel&& kernel,
            std::index_sequence<I...>
        ) const {
            using K = std::remove_cvref_t<Kernel>;

            std::forward<Kernel>(kernel)(
                container_ptr->template at_id<
                    K::Read,
                    K::Write
                >(ids[I])...
            );
        }
    };


    // @Brief Builds phases of independent/non-conflicting Topological Batches (collection of id-id interactions)
    template <typename Node = ParticleID>
    std::vector<std::vector<utility::graph::EdgeList<Node>>> build_concurrent_phases(
        const std::vector<utility::graph::EdgeList<Node>> & input_topologies,
        const size_t max_partition_size,
        const size_t min_batches_threshold
    ) {
        using Pairs = utility::graph::EdgeList<Node>;
        using Set = std::vector<Node>;

        // find all connected components (atomics) for all input topologies
        // pairs of particles inside an atomic batch must be scheduled together
        std::vector<Pairs> atomics;
        for (const auto& edges : input_topologies) {
            if (edges.empty()) continue;

            const auto components = utility::graph::find_connected_components(edges);
            for (auto& comp : components) {
                atomics.push_back(comp);
            }
        }

        const size_t num_batches = atomics.size();
        if (num_batches == 0) return {};

        // Build a conflict graph (intersection graph i.e. which atomics have common nodes)
        // First get all nodes in each atomic
        std::vector<Set> node_sets;
        node_sets.reserve(num_batches);
        for (const auto& atomic_batch : atomics) {
            node_sets.push_back(utility::graph::get_unique_nodes(atomic_batch));
        }
        // then build conflict graph
        const auto conflict_graph = utility::graph::build_intersection_graph(node_sets);

        // build phases (sets) of independent atomics (i.e. concurrently schedulable)
        // by using a greedy partitioning (coloring) algorithm on the conflict graph
        // We first sort the nodes of the conflict graph (the atomics) by size
        std::vector<size_t> processing_order(num_batches);
        std::iota(processing_order.begin(), processing_order.end(), 0);
        std::ranges::sort(processing_order, [&](const size_t a, const size_t b) {
            const size_t degree_a = conflict_graph[a].size();
            const size_t degree_b = conflict_graph[b].size();
            if (degree_a != degree_b) return degree_a > degree_b;
            return atomics[a].size() > atomics[b].size();
        });

        // the result of the partitioning algorithm is a list of sets of indices. The indices index into the atomics vector
        const auto packed_phases = utility::graph::greedy_independent_partitions(
            conflict_graph, processing_order, max_partition_size
        );

        // build phases with particle pairs back from "index" phases (packed_phases)
        std::vector<std::vector<Pairs>> final_phases;
        std::vector<size_t> sequential_fallback;

        for (const auto& phase_indices : packed_phases) {
            // a phase must have a certain number of concurrently executable batches
            // else we execute them sequentially
            if (phase_indices.size() < min_batches_threshold) {
                sequential_fallback.insert(
                    sequential_fallback.end(),
                    phase_indices.begin(),
                    phase_indices.end()
                );
                continue;
            }

            std::vector<Pairs> current_phase;
            current_phase.reserve(phase_indices.size());
            for (const size_t idx : phase_indices) {
                current_phase.push_back(std::move(atomics[idx]));
            }
            final_phases.push_back(std::move(current_phase));
        }

        // Add the sequential fallback phase
        if (!sequential_fallback.empty()) {
            Pairs seq_batch;

            // pre allocate space
            size_t total_pairs = 0;
            for (const size_t idx : sequential_fallback) {
                total_pairs += atomics[idx].size();
            }
            seq_batch.reserve(total_pairs);

            // Flatten all under-utilized atomics into seq_batch
            for (const size_t idx : sequential_fallback) {
                seq_batch.insert(
                    seq_batch.end(),
                    std::make_move_iterator(atomics[idx].begin()),
                    std::make_move_iterator(atomics[idx].end())
                );
            }

            std::vector<Pairs> seq_phase {std::move(seq_batch)};
            final_phases.push_back(std::move(seq_phase));
        }

        return final_phases;
    }
}
