#include <gtest/gtest.h>
#include <algorithm>

#include "april/utility/graph.hpp"
using namespace april::utility::graph;

// Helper to sort a vector of vectors for deterministic comparisons
template <typename T>
void normalize_2d_vector(std::vector<std::vector<T>>& vec) {
    for (auto& inner : vec) {
        std::ranges::sort(inner);
    }
    std::ranges::sort(vec, [](const auto& a, const auto& b) {
        if (a.empty() || b.empty()) return a.size() < b.size();
        return a[0] < b[0];
    });
}

// Helper to sort components (EdgeLists)
template <typename Node>
void normalize_components(std::vector<EdgeList<Node>>& components) {
    for (auto& comp : components) {
        std::ranges::sort(comp);
    }
    std::ranges::sort(components, [](const auto& a, const auto& b) {
        if (a.empty() || b.empty()) return a.size() < b.size();
        return a[0] < b[0]; // Compare by first pair
    });
}

TEST(GraphUtilityTest, FindConnectedComponents) {
    EdgeList<int> edges = {
        {1, 2}, {2, 3}, // Component 1
        {4, 5}, {5, 6}, // Component 2
        {7, 8}          // Component 3
    };

    auto components = find_connected_components(edges);
    normalize_components(components);

    ASSERT_EQ(components.size(), 3);

    const EdgeList<int> expected_comp1 = {{1, 2}, {2, 3}};
    const EdgeList<int> expected_comp2 = {{4, 5}, {5, 6}};
    const EdgeList<int> expected_comp3 = {{7, 8}};

    EXPECT_EQ(components[0], expected_comp1);
    EXPECT_EQ(components[1], expected_comp2);
    EXPECT_EQ(components[2], expected_comp3);
}

TEST(GraphUtilityTest, FindConnectedComponentsEmpty) {
    EdgeList<int> edges = {};
    auto components = find_connected_components(edges);
    EXPECT_TRUE(components.empty());
}

TEST(GraphUtilityTest, GetUniqueNodes) {
    EdgeList<int> edges = {{1, 2}, {2, 3}, {1, 4}, {5, 5}};
    auto nodes = get_unique_nodes(edges);

    // get_unique_nodes inherently sorts the output
    std::vector<int> expected = {1, 2, 3, 4, 5};
    EXPECT_EQ(nodes, expected);
}

TEST(GraphUtilityTest, BuildIntersectionGraph) {
    std::vector<std::vector<int>> node_sets = {
        {1, 2, 3}, // Set 0
        {3, 4},    // Set 1 (intersects 0 at 3)
        {5, 6},    // Set 2 (disjoint)
        {4, 6}     // Set 3 (intersects 1 at 4, intersects 2 at 6)
    };

    auto adj_list = build_intersection_graph(node_sets);

    ASSERT_EQ(adj_list.size(), 4);

    // Internal adjacency lists are sorted in your implementation
    EXPECT_EQ(adj_list[0], std::vector<size_t>{1});
    EXPECT_EQ(adj_list[1], (std::vector<size_t>{0, 3}));
    EXPECT_EQ(adj_list[2], std::vector<size_t>{3});
    EXPECT_EQ(adj_list[3], (std::vector<size_t>{1, 2}));
}

TEST(GraphUtilityTest, GreedyIndependentPartitions) {
    // A square graph: 0-1-3-2-0
    AdjacencyList<size_t> adj_list = {
        {1, 2}, // Node 0 connected to 1, 2
        {0, 3}, // Node 1 connected to 0, 3
        {0, 3}, // Node 2 connected to 0, 3
        {1, 2}  // Node 3 connected to 1, 2
    };

    std::vector<size_t> order = {0, 1, 2, 3};
    auto partitions = greedy_independent_partitions(adj_list, order);

    normalize_2d_vector(partitions);

    ASSERT_EQ(partitions.size(), 2);
    // Nodes 0 and 3 are independent. Nodes 1 and 2 are independent.
    EXPECT_EQ(partitions[0], (std::vector<size_t>{0, 3}));
    EXPECT_EQ(partitions[1], (std::vector<size_t>{1, 2}));
}

TEST(GraphUtilityTest, GreedyIndependentPartitionsWithMaxSize) {
    // Star graph: Node 0 is center, connected to 1, 2, 3
    AdjacencyList<size_t> adj_list = {
        {1, 2, 3}, // Node 0
        {0},       // Node 1
        {0},       // Node 2
        {0}        // Node 3
    };

    std::vector<size_t> order = {0, 1, 2, 3};
    // Force a max partition size of 2
    auto partitions = greedy_independent_partitions(adj_list, order, 2);

    normalize_2d_vector(partitions);

    // Node 0 is alone. Nodes 1,2,3 want to be together, but max size is 2.
    // So we expect 3 partitions: {0}, {1, 2}, {3}
    ASSERT_EQ(partitions.size(), 3);
    EXPECT_EQ(partitions[0], std::vector<size_t>{0});
    EXPECT_EQ(partitions[1], (std::vector<size_t>{1, 2}));
    EXPECT_EQ(partitions[2], std::vector<size_t>{3});
}


TEST(GraphUtilityTest, GetUniqueNodesEmpty) {
    EdgeList<int> edges = {};
    auto nodes = get_unique_nodes(edges);
    EXPECT_TRUE(nodes.empty());
}

TEST(GraphUtilityTest, BuildIntersectionGraphEmpty) {
    std::vector<std::vector<int>> node_sets = {};
    auto adj_list = build_intersection_graph(node_sets);
    EXPECT_TRUE(adj_list.empty());
}

TEST(GraphUtilityTest, GreedyIndependentPartitionsEmpty) {
    AdjacencyList<size_t> adj_list = {};
    std::vector<size_t> order = {};
    auto partitions = greedy_independent_partitions(adj_list, order);
    EXPECT_TRUE(partitions.empty());
}

TEST(GraphUtilityTest, BuildIntersectionGraphMultipleOverlaps) {
    const std::vector<std::vector<int>> node_sets = {
        {1, 2, 3}, // Set 0
        {2, 3, 4}  // Set 1 (Shares both 2 and 3 with Set 0)
    };

    auto adj_list = build_intersection_graph(node_sets);

    ASSERT_EQ(adj_list.size(), 2);
    // Even though they share two nodes, there should only be one edge between them
    EXPECT_EQ(adj_list[0], std::vector<size_t>{1});
    EXPECT_EQ(adj_list[1], std::vector<size_t>{0});
}

TEST(GraphUtilityTest, GreedyIndependentPartitionsNoEdges) {
    // 4 nodes, no edges between them
    AdjacencyList<size_t> adj_list(4);

    std::vector<size_t> order = {0, 1, 2, 3};
    auto partitions = greedy_independent_partitions(adj_list, order, 2);

    normalize_2d_vector(partitions);

    // Because max size is 2 and no nodes are forbidden, we expect exactly 2 partitions of size 2.
    ASSERT_EQ(partitions.size(), 2);
    EXPECT_EQ(partitions[0], (std::vector<size_t>{0, 1}));
    EXPECT_EQ(partitions[1], (std::vector<size_t>{2, 3}));
}