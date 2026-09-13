#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <variant>

#include "april/base/types.hpp"
#include "april/interactions/interaction.hpp"
#include "april/interactions/coulomb.hpp"
#include "april/interactions/interaction_table.hpp"
#include "constant_force.h"

using namespace april;


// 1. Define the InteractionVariant and InteractionTable types for the test
// Must include InteractionSentinel and NoInteraction as per internal requirements
using TestInteractionVariant = std::variant<interaction::internal::InteractionSentinel, ConstantForce, NoInteraction>;
using InteractionTable = interaction::internal::InteractionTable<TestInteractionVariant>;
using TypeInfo = interaction::internal::TypeInteraction<TestInteractionVariant>;
using IdInfo = interaction::internal::IdInteraction<TestInteractionVariant>;


TEST(InteractionManagerTest, EmptyBuild) {
    EXPECT_NO_THROW(InteractionTable({}, {}, {}, {}));

    const InteractionTable interaction_table({}, {}, {}, {});

    auto schema = interaction_table.generate_interaction_map();
    EXPECT_TRUE(schema.interactions.empty());
}

TEST(InteractionManagerTest, MaxCutoffCalculation) {
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(1, 1, 1, 1.5));
    info.emplace_back(1, 1, ConstantForce(2, 2, 2, 2.5));

    const std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};

    const InteractionTable interaction_table(info, {}, type_map, {});

    // Verify via Schema
    auto schema = interaction_table.generate_interaction_map();

    double max_cut = 0;
    for(auto& p : schema.interactions) {
        if(p.is_active) max_cut = std::max(max_cut, p.cutoff);
    }

    EXPECT_DOUBLE_EQ(max_cut, 2.5);
}

TEST(InteractionManagerTest, TypeBasedLookup) {
    std::vector<TypeInfo> info;
    // Define interactions:
    // 0-0: (4,5,6)
    // 1-1: (1,2,3)
    // 0-1: (7,8,9)
    info.emplace_back(0, 0, ConstantForce(4, 5, 6));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3));
    info.emplace_back(0, 1, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    InteractionTable interaction_table(info, {}, type_map, {});

    // Helper to run dispatch and extract the force vector
    auto eval_type = [&](ParticleType t1, ParticleType t2) {
        vec3 result{0,0,0};
        interaction_table.dispatch(t1, t2, [&](const auto& interaction) {
            if constexpr (std::is_same_v<std::decay_t<decltype(interaction)>, ConstantForce>) {
                 result = interaction.v; // UPDATED: Using .v
            }
        });
        return result;
    };

    EXPECT_EQ(eval_type(0, 0), vec3(4, 5, 6));
    EXPECT_EQ(eval_type(1, 1), vec3(1, 2, 3));

    // Check symmetry
    EXPECT_EQ(eval_type(0, 1), vec3(7, 8, 9));
    EXPECT_EQ(eval_type(1, 0), vec3(7, 8, 9));
}

TEST(InteractionManagerTest, IdBasedLookup) {
    std::vector<IdInfo> id_info;
    std::vector<TypeInfo> type_info;

    // Type interaction (0,0) is Zero
    type_info.emplace_back(0, 0, ConstantForce(0, 0, 0));

    // ID interaction (42, 99) is (7,8,9)
    // internal map: 42->0, 99->1
    id_info.emplace_back(42, 99, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}};
    std::unordered_map<ParticleID, ParticleID> id_map{{42, 0}, {99, 1}, {100, 2}};
    const InteractionTable interaction_table(type_info, id_info, type_map, id_map);

    // 1. Check ID interaction existence
    EXPECT_TRUE(interaction_table.has_id_interaction(0, 1));
    EXPECT_TRUE(interaction_table.has_id_interaction(1, 0));
    EXPECT_TRUE(interaction_table.has_id_interaction(0, 0));
    EXPECT_TRUE(interaction_table.has_id_interaction(1, 1));
    EXPECT_FALSE(interaction_table.has_id_interaction(0, 2));
    EXPECT_FALSE(interaction_table.has_id_interaction(2, 2));

    // 2. Dispatch ID
    auto eval_id = [&](ParticleID id1, ParticleID id2) {
        vec3 result{0,0,0};
        interaction_table.dispatch_id(id1, id2, [&](const auto& interaction) {
            if constexpr (std::is_same_v<std::decay_t<decltype(interaction)>, ConstantForce>) {
                 result = interaction.v; // UPDATED: Using .v
            }
        });
        return result;
    };

    // Note: Inputs to dispatch_id are Implementation IDs
    EXPECT_EQ(eval_id(0, 1), vec3(7, 8, 9));
    EXPECT_EQ(eval_id(1, 0), vec3(7, 8, 9));

    // 3. Verify Schema Topology
    auto schema = interaction_table.generate_interaction_map();

    // We expect the schema to have recorded the ID usage for this pair
    bool found_id_link = false;
    for(const auto& prop : schema.interactions) {
        for(const auto& pair : prop.used_by_ids) {
            // Check for pair (0, 1) or (1, 0)
            if ((pair.first == 0 && pair.second == 1) || (pair.first == 1 && pair.second == 0)) {
                found_id_link = true;
            }
        }
    }
    EXPECT_TRUE(found_id_link) << "Schema should record the ID usage for (0,1)";
}

TEST(InteractionManagerTest, MixingInteractions) {
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(4, 5, 6));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3));
    // Missing (0,1) -> Should trigger mixing!

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    InteractionTable interaction_table(info, {}, type_map, {});

    auto eval_type = [&](ParticleType t1, ParticleType t2) {
        vec3 result{0,0,0};
        interaction_table.dispatch(t1, t2, [&](const auto& interaction) {
            if constexpr (std::is_same_v<std::decay_t<decltype(interaction)>, ConstantForce>) {
                 result = interaction.v;
            }
        });
        return result;
    };

    EXPECT_EQ(eval_type(0, 0), vec3(4, 5, 6));
    EXPECT_EQ(eval_type(1, 1), vec3(1, 2, 3));

    // UPDATED: ConstantForce::mix sums the vectors
    // (4,5,6) + (1,2,3) = (5, 7, 9)
    vec3 expected(5.0, 7.0, 9.0);
    EXPECT_EQ(eval_type(0, 1), expected);
    EXPECT_EQ(eval_type(1, 0), expected);
}

TEST(InteractionManagerTest, SchemaDeduplication) {
    std::vector<TypeInfo> info;
    // Two different pairs use identical interactions
    info.emplace_back(0, 0, ConstantForce(1, 0, 0));
    info.emplace_back(1, 1, ConstantForce(1, 0, 0)); // Same as 0-0
    info.emplace_back(0, 1, ConstantForce(2, 0, 0)); // Different

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    InteractionTable interaction_table(info, {}, type_map, {});

    auto schema = interaction_table.generate_interaction_map();

    // We expect exactly 2 unique interactions in the palette:
    // 1. ConstantForce(1,0,0) [used by 0-0 and 1-1]
    // 2. ConstantForce(2,0,0) [used by 0-1]
    EXPECT_EQ(schema.interactions.size(), 2);
}











