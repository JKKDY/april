#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <variant>

#include "april/common.hpp"
#include "april/forces/force.hpp"
#include "april/forces/coulomb.hpp"
#include "april/forces/force_table.hpp"
#include "constant_force.h" // Assumes your provided struct is here

using namespace april;
using namespace april::force;

// 1. Define the ForceVariant and ForceTable types for the test
// Must include ForceSentinel and NoForce as per internal requirements
using TestForceVariant = std::variant<internal::ForceSentinel, ConstantForce, NoForce>;
using ForceTable = internal::ForceTable<TestForceVariant>;
using TypeInfo = internal::TypeInteraction<TestForceVariant>;
using IdInfo = internal::IdInteraction<TestForceVariant>;


TEST(InteractionManagerTest, EmptyBuild) {
    EXPECT_NO_THROW(ForceTable({}, {}, {}, {}));

    const ForceTable force_table({}, {}, {}, {});

    auto schema = force_table.generate_schema();
    EXPECT_TRUE(schema.interactions.empty());
}

TEST(InteractionManagerTest, MaxCutoffCalculation) {
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(1, 1, 1, 1.5));
    info.emplace_back(1, 1, ConstantForce(2, 2, 2, 2.5));

    const std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};

    const ForceTable force_table(info, {}, type_map, {});

    // Verify via Schema
    auto schema = force_table.generate_schema();

    double max_cut = 0;
    for(auto& p : schema.interactions) {
        if(p.is_active) max_cut = std::max(max_cut, p.cutoff);
    }

    EXPECT_DOUBLE_EQ(max_cut, 2.5);
}

TEST(InteractionManagerTest, TypeBasedLookup) {
    std::vector<TypeInfo> info;
    // Define forces:
    // 0-0: (4,5,6)
    // 1-1: (1,2,3)
    // 0-1: (7,8,9)
    info.emplace_back(0, 0, ConstantForce(4, 5, 6));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3));
    info.emplace_back(0, 1, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    ForceTable force_table(info, {}, type_map, {});

    // Helper to run dispatch and extract the force vector
    auto eval_type = [&](ParticleType t1, ParticleType t2) {
        vec3 result{0,0,0};
        force_table.dispatch(t1, t2, [&](const auto& force) {
            if constexpr (std::is_same_v<std::decay_t<decltype(force)>, ConstantForce>) {
                 result = force.v; // UPDATED: Using .v
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

    // Type Force (0,0) is Zero
    type_info.emplace_back(0, 0, ConstantForce(0, 0, 0));

    // ID Force (42, 99) is (7,8,9)
    // internal map: 42->0, 99->1
    id_info.emplace_back(42, 99, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}};
    std::unordered_map<ParticleID, ParticleID> id_map{{42, 0}, {99, 1}, {100, 2}};
    const ForceTable force_table(type_info, id_info, type_map, id_map);

    // 1. Check ID interaction existence
    EXPECT_TRUE(force_table.has_id_force(0, 1));
    EXPECT_TRUE(force_table.has_id_force(1, 0));
    EXPECT_TRUE(force_table.has_id_force(0, 0));
    EXPECT_TRUE(force_table.has_id_force(1, 1));
    EXPECT_FALSE(force_table.has_id_force(0, 2));
    EXPECT_FALSE(force_table.has_id_force(2, 2));

    // 2. Dispatch ID
    auto eval_id = [&](ParticleID id1, ParticleID id2) {
        vec3 result{0,0,0};
        force_table.dispatch_id(id1, id2, [&](const auto& force) {
            if constexpr (std::is_same_v<std::decay_t<decltype(force)>, ConstantForce>) {
                 result = force.v; // UPDATED: Using .v
            }
        });
        return result;
    };

    // Note: Inputs to dispatch_id are Implementation IDs
    EXPECT_EQ(eval_id(0, 1), vec3(7, 8, 9));
    EXPECT_EQ(eval_id(1, 0), vec3(7, 8, 9));

    // 3. Verify Schema Topology
    auto schema = force_table.generate_schema();

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

TEST(InteractionManagerTest, MixingForces) {
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(4, 5, 6));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3));
    // Missing (0,1) -> Should trigger mixing!

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    ForceTable force_table(info, {}, type_map, {});

    auto eval_type = [&](ParticleType t1, ParticleType t2) {
        vec3 result{0,0,0};
        force_table.dispatch(t1, t2, [&](const auto& force) {
            if constexpr (std::is_same_v<std::decay_t<decltype(force)>, ConstantForce>) {
                 result = force.v; // UPDATED: Using .v
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
    // Two different pairs use Identical Forces
    info.emplace_back(0, 0, ConstantForce(1, 0, 0));
    info.emplace_back(1, 1, ConstantForce(1, 0, 0)); // Same as 0-0
    info.emplace_back(0, 1, ConstantForce(2, 0, 0)); // Different

    std::unordered_map<ParticleType, ParticleType> type_map{{0, 0}, {1, 1}};
    ForceTable force_table(info, {}, type_map, {});

    auto schema = force_table.generate_schema();

    // We expect exactly 2 unique interactions in the palette:
    // 1. Force(1,0,0) [used by 0-0 and 1-1]
    // 2. Force(2,0,0) [used by 0-1]
    EXPECT_EQ(schema.interactions.size(), 2);
}