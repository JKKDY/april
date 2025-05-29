// interaction_manager_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <memory>

#include "april/env/interaction.h"
#include "april/env/particle.h"
#include "april/common.h"

using namespace april;
using namespace april::env;
using InteractionManager = impl::InteractionManager;
using InteractionInfo = impl::InteractionInfo;

// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final : Force {
    vec3 v;
    ConstantForce(double x, double y, double z, double cutoff = -1) : v{x,y,z} {
        cutoff_radius = cutoff;
    }
    vec3 operator()(const impl::Particle&, const impl::Particle&, const vec3&) const noexcept override {
        return v;
    }
    std::unique_ptr<Force> mix(const Force* other) const override {
        auto* o = dynamic_cast<const ConstantForce*>(other);
        if (!o) throw std::invalid_argument("mix mismatch");
        return std::make_unique<ConstantForce>(
            v.x + o->v.x,
            v.y + o->v.y,
            v.z + o->v.z,
            std::max(cutoff_radius, o->cutoff_radius)
        );
    }
};

// Helper to make a dummy particle
static impl::Particle make_particle(impl::ParticleType type, impl::ParticleID id, double mass=1.0, vec3 pos={0,0,0}) {
    return {
    /* index       */ id,
    /* id          */ id,
    /* position    */ pos,
    /* velocity    */ vec3{0,0,0},
    /* mass        */ mass,
    /* type        */ type,
    /* state       */ impl::Particle::State::ALIVE
    };
}

TEST(InteractionManagerTest, EmptyBuild) {
    InteractionManager mgr;
    std::vector<InteractionInfo> info;
    // empty maps fine
    EXPECT_NO_THROW(mgr.build(info, {}, {}));
    EXPECT_EQ(mgr.get_max_cutoff(), 0.0);
}

TEST(InteractionManagerTest, MaxCutoffCalculation) {
    InteractionManager mgr;

    // two user type‐based interactions with cutoffs 1.5 and 2.5
    std::vector<InteractionInfo> info;
    info.emplace_back(true, std::pair{0,0},std::make_unique<ConstantForce>(1,1,1,1.5));
    info.emplace_back(true, std::pair{1,1},std::make_unique<ConstantForce>(2,2,2,2.5));

    std::unordered_map<ParticleType, impl::ParticleType> type_map {{0,0},{1,1}};

    EXPECT_NO_THROW(mgr.build(info, type_map, {}));
    EXPECT_DOUBLE_EQ(mgr.get_max_cutoff(), 2.5);
}

TEST(InteractionManagerTest, TypeBasedLookup) {
    InteractionManager mgr;

    std::vector<InteractionInfo> info;
    info.emplace_back(true, std::pair{0,0},std::make_unique<ConstantForce>(4,5,6, -1));
    info.emplace_back(true, std::pair{1,1},std::make_unique<ConstantForce>(1,2,3, -1));
    info.emplace_back(true, std::pair{0,1},std::make_unique<ConstantForce>(7,8,9, -1));

    std::unordered_map<ParticleType, impl::ParticleType> type_map {{0,0},{1,1}};
    mgr.build(info, type_map, {});

    auto p0 = make_particle(0, 10, 1.0, {0,0,0});
    auto p1 = make_particle(1, 11, 1.0, {1,1,1});

    vec3 f1 = mgr.evaluate(p0, p0); // in actual use a particle should never interact with itself ofc
    EXPECT_EQ(f1, vec3(4,5,6));

    vec3 f2 = mgr.evaluate(p1, p1);
    EXPECT_EQ(f2, vec3(1,2,3));

    vec3 f3 = mgr.evaluate(p0, p1);
    vec3 f4 = mgr.evaluate(p1, p0);
    EXPECT_EQ(f3, f4);
    EXPECT_EQ(f3, vec3(7,8,9));
}

TEST(InteractionManagerTest, IdBasedLookup) {
    InteractionManager mgr;
    // one id‐based entry for (42,99)
    std::vector<InteractionInfo> info;
    info.emplace_back(false, std::pair{42,99}, std::make_unique<ConstantForce>(7,8,9));
    info.emplace_back(true, std::pair{0,0},std::make_unique<NoForce>());


    std::unordered_map<ParticleType, impl::ParticleType> type_map { {0,0}};
    std::unordered_map<ParticleID, impl::ParticleID> id_map {{42,0},{99,1}};
    mgr.build(info, type_map, id_map);

    auto p1 = make_particle(0, 0);
    auto p2 = make_particle(0, 1);
    auto p3 = make_particle(0, 2);

    vec3 f1 = mgr.evaluate(p1, p2);
    vec3 f2 = mgr.evaluate(p2, p1);
    EXPECT_EQ(f1, f2);
    EXPECT_EQ(f1, vec3(7,8,9));

    vec3 f = mgr.evaluate(p3, p3);

    EXPECT_EQ(f, vec3{});
}


TEST(InteractionManagerTest, MixingForces) {
    InteractionManager mgr;

    std::vector<InteractionInfo> info;
    info.emplace_back(true, std::pair{0,0},std::make_unique<ConstantForce>(4,5,6, -1));
    info.emplace_back(true, std::pair{1,1},std::make_unique<ConstantForce>(1,2,3, -1));

    std::unordered_map<ParticleType, impl::ParticleType> type_map {{0,0},{1,1}};
    mgr.build(info, type_map, {});

    auto p0 = make_particle(0, 10, 1.0, {0,0,0});
    auto p1 = make_particle(1, 11, 1.0, {1,1,1});

    vec3 f1 = mgr.evaluate(p0, p0);
    EXPECT_EQ(f1, vec3(4,5,6));

    vec3 f2 = mgr.evaluate(p1, p1);
    EXPECT_EQ(f2, vec3(1,2,3));

    vec3 f3 = mgr.evaluate(p0, p1);
    vec3 f4 = mgr.evaluate(p1, p0);
    EXPECT_EQ(f3, f4);
    EXPECT_EQ(f3, vec3(5,7,9));
}
