// interaction_manager_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <memory>

# include "april/april.h"
#include "constant_force.h"

using namespace april;


// Helper to make a dummy particle
static env::internal::Particle make_particle(env::internal::ParticleType type, env::internal::ParticleID id, double mass=1.0, vec3 pos={0,0,0}) {
    return {
    /* id          */ id,
    /* position    */ pos,
    /* velocity    */ vec3{0,0,0},
    /* mass        */ mass,
    /* type        */ type,
    /* state       */ env::internal::Particle::State::ALIVE
    };
}


// Use an environment that supports ConstantForce
using Env = Environment<force::ForcePack<ConstantForce>, boundary::BoundaryPack<>, controller::ControllerPack<>, field::FieldPack<>>;
using FT  = Env::traits::force_table_t;
using TypeInfo = force::internal::TypeInteraction<Env::traits::force_variant_t>; // variant<ConstantForce>
using IdInfo = force::internal::IdInteraction<Env::traits::force_variant_t>; // variant<ConstantForce>


TEST(InteractionManagerTest, EmptyBuild) {
    // empty maps are fine
    EXPECT_NO_THROW(FT({}, {}, {}, {}));

    const FT mgr({}, {}, {}, {});
    EXPECT_DOUBLE_EQ(mgr.get_max_cutoff(), 0.0);
}

TEST(InteractionManagerTest, MaxCutoffCalculation) {

    // two type-based interactions with cutoffs 1.5 and 2.5
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(1, 1, 1, 1.5));
    info.emplace_back(1, 1, ConstantForce(2, 2, 2, 2.5));

    const std::unordered_map<ParticleType, env::internal::ParticleType> type_map{{0, 0}, {1, 1}};

    const FT mgr(info, {}, type_map, {});
    EXPECT_NO_THROW(FT(info, {}, type_map, {}));
    EXPECT_DOUBLE_EQ(mgr.get_max_cutoff(), 2.5);
}

TEST(InteractionManagerTest, TypeBasedLookup) {
    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(4, 5, 6, -1));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3, -1));
    info.emplace_back(0, 1, ConstantForce(7, 8, 9, -1));

    std::unordered_map<ParticleType, env::internal::ParticleType> type_map{{0, 0}, {1, 1}};
    FT mgr(info, {}, type_map, {});

    auto p0 = make_particle(0, 10, 1.0, {0, 0, 0});
    auto p1 = make_particle(1, 11, 1.0, {1, 1, 1});

    vec3 f1 = mgr.evaluate(p0, p0); // self-interaction (just for lookup test)
    EXPECT_EQ(f1, vec3(4, 5, 6));

    vec3 f2 = mgr.evaluate(p1, p1);
    EXPECT_EQ(f2, vec3(1, 2, 3));

    vec3 f3 = mgr.evaluate(p0, p1);
    vec3 f4 = mgr.evaluate(p1, p0);
    EXPECT_EQ(f3, f4);
    EXPECT_EQ(f3, vec3(7, 8, 9));
}

TEST(InteractionManagerTest, IdBasedLookup) {
    std::vector<IdInfo> id_info;
    std::vector<TypeInfo> type_info;

    // Provide a zero type-force for (0,0) so evaluate never hits NullForce
    type_info.emplace_back(0, 0, ConstantForce(0, 0, 0));
    // one id-based entry for (42,99)
    id_info.emplace_back(42, 99, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, env::internal::ParticleType> type_map{{0, 0}};
    std::unordered_map<ParticleID, env::internal::ParticleID> id_map{{42, 0}, {99, 1}};
    const FT mgr (type_info, id_info, type_map, id_map);

    auto p1 = make_particle(0, 0);
    auto p2 = make_particle(0, 1);
    auto p3 = make_particle(0, 2);

    vec3 f1 = mgr.evaluate(p1, p2);
    vec3 f2 = mgr.evaluate(p2, p1);
    EXPECT_EQ(f1, f2);
    EXPECT_EQ(f1, vec3(7, 8, 9));

    // no id interaction for (2,2), and type force is zero; expect zero
    vec3 f = mgr.evaluate(p3, p3);
    EXPECT_EQ(f, vec3{});
}


TEST(InteractionManagerTest, MixingForces) {

    std::vector<TypeInfo> info;
    info.emplace_back(0, 0, ConstantForce(4, 5, 6, -1));
    info.emplace_back(1, 1, ConstantForce(1, 2, 3, -1));

    std::unordered_map<ParticleType, env::internal::ParticleType> type_map{{0, 0}, {1, 1}};
    FT mgr (info, {}, type_map, {});

    auto p0 = make_particle(0, 10, 1.0, {0, 0, 0});
    auto p1 = make_particle(1, 11, 1.0, {1, 1, 1});

    vec3 f1 = mgr.evaluate(p0, p0);
    EXPECT_EQ(f1, vec3(4, 5, 6));

    vec3 f2 = mgr.evaluate(p1, p1);
    EXPECT_EQ(f2, vec3(1, 2, 3));

    vec3 f3 = mgr.evaluate(p0, p1);
    vec3 f4 = mgr.evaluate(p1, p0);
    EXPECT_EQ(f3, f4);
    EXPECT_EQ(f3, vec3(5, 7, 9)); // assuming your container mixes type (0,0) with (1,1) for cross-type
}
