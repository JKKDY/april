#pragma once
#include <ranges>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>


namespace april::utility::graph {

    // Edge list representation of a graph
    template <typename Node>
    using EdgeList = std::vector<std::pair<Node, Node>>;

    // Index-based adjacency list representation of a graph
    template <typename Node = size_t>
    using AdjacencyList = std::vector<std::vector<Node>>;


    // @brief Find all connected components in an edge list graph using union-find
    template <typename Node>
    std::vector<EdgeList<Node>> find_connected_components(const EdgeList<Node>& edges) {
        if (edges.empty()) return {};

        struct {
            std::unordered_map<Node, Node> parent;

            // find representative of the connected component containing id (i.e. the root parent of id)
            Node find(const Node i) {
                if (!parent.contains(i)) { // no parent? set i to be its own parent
                    parent[i] = i;
                    return i;
                }
                if (parent[i] == i) {
                    return i;
                }
                // path compression for shorter look-ups
                return parent[i] = find(parent[i]);
            }

            // unite two connected components containing i and j
            void unite(const Node i, const Node j) {
                const Node root_i = find(i);
                const Node root_j = find(j);
                if (root_i != root_j) {
                    parent[root_i] = root_j;
                }
            }
        } dsu; // disjoint-set-union: standard Union-Find data structure to identify connected components

        // find all connected components
        for (const auto& pair : edges) {
            dsu.unite(pair.first, pair.second);
        }

        // group the pairs by their component's representative ID
        std::unordered_map<Node, EdgeList<Node>> component_map;
        for (const auto& pair : edges) {
            // Both first and second share the same root at this point
            Node representative = dsu.find(pair.first);
            component_map[representative].push_back(pair);
        }

        // Move grouped batches into a vector
        std::vector<EdgeList<Node>> connected_components;
        connected_components.reserve(component_map.size());
        for (auto& batch : component_map | std::views::values) {
            connected_components.push_back(std::move(batch));
        }

        return connected_components;
    }


    // @brief Return all (unique) nodes given an edge-list graph
    template <typename Node>
    std::vector<Node> get_unique_nodes(const EdgeList<Node>& edges) {
        std::vector<Node> nodes;
        nodes.reserve(edges.size() * 2);

        for (const auto& [u, v] : edges) {
            nodes.push_back(u);
            nodes.push_back(v);
        }

        std::ranges::sort(nodes);
        auto [first, last] = std::ranges::unique(nodes);
        nodes.erase(first, last);

        return nodes;
    }


    // @brief Build an Intersection Graph (Adjacency List) from collections of nodes
    template <typename Node>
    AdjacencyList<size_t> build_intersection_graph(const std::vector<std::vector<Node>>& node_sets) {
        const size_t num_sets = node_sets.size();
        AdjacencyList<size_t> adj_list(num_sets);

        // Map each node to the indices of the sets it belongs to (inverted index)
        std::unordered_map<Node, std::vector<size_t>> inverted_index;
        for (size_t i = 0; i < num_sets; ++i) {
            for (const Node& n : node_sets[i]) {
                inverted_index[n].push_back(i);
            }
        }

        // Build adjacency list: sets sharing a node are adjacent
        for (const auto& shared_indices : inverted_index | std::views::values) {
            for (size_t i = 0; i < shared_indices.size(); ++i) {
                for (size_t j = i + 1; j < shared_indices.size(); ++j) {
                    size_t c1 = shared_indices[i];
                    size_t c2 = shared_indices[j];

                    adj_list[c1].push_back(c2);
                    adj_list[c2].push_back(c1);
                }
            }
        }

        // Clean up redundant edges caused by sets sharing multiple nodes
        for (auto& neighbors : adj_list) {
            std::ranges::sort(neighbors);
            auto [first, last] = std::ranges::unique(neighbors);
            neighbors.erase(first, last);
        }

        return adj_list;
    }


    // @brief Partitions vertices into independent sets/colors (no edge between vertices of a given color) given an AdjacencyList.
    // Respects a given processing order and an optional maximum partition size.
    template <typename Node = size_t>
    std::vector<std::vector<Node>> greedy_independent_partitions(
         const AdjacencyList<Node>& adj_list,
         const std::vector<Node>& processing_order,
         size_t max_partition_size = std::numeric_limits<size_t>::max())
    {
        struct Partition {
            std::vector<Node> nodes;
            std::unordered_set<Node> forbidden; // Neighbors of nodes in this partition
        };

        std::vector<Partition> partitions;

        for (Node u : processing_order) {
            bool placed = false;

            for (auto& partition : partitions) {
                // if partition already large enough skip and try to add to a different partition (or make a new one)
                if (partition.nodes.size() >= max_partition_size) continue;

                // if u is not in the connected to the partition add it and add its neighbors to the forbidden list
                if (!partition.forbidden.contains(u)) {
                    partition.nodes.push_back(u);
                    for (Node neighbor : adj_list[u]) {
                        partition.forbidden.insert(neighbor);
                    }
                    placed = true;
                    break;
                }
            }

            // create a new partition if u doesn't fit in any of the other partitions
            if (!placed) {
                Partition new_part;
                new_part.nodes.push_back(u);
                for (Node neighbor : adj_list[u]) {
                    new_part.forbidden.insert(neighbor);
                }
                partitions.push_back(std::move(new_part));
            }
        }

        // convert std::vector<Partition> to std::vector<std::vector<Node>> and return result
        std::vector<std::vector<Node>> result;
        result.reserve(partitions.size());
        for (auto& part : partitions) {
            result.push_back(std::move(part.nodes));
        }

        return result;
    }
}