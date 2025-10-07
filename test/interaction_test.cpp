// interaction_manager_test.cpp
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <memory>

# include "april/april.h"
using namespace april;

template <class Env>
using InteractionManager = env::impl::InteractionManager<Env>;

// A tiny force that returns a constant vector and mixes by summing
struct ConstantForce final {
    vec3 v;
    double cutoff_radius;
    ConstantForce(double x, double y, double z, double cutoff = -1) : v{x,y,z} {
        cutoff_radius = cutoff;
    }
    vec3 operator()(const env::impl::Particle&, const env::impl::Particle&, const vec3&) const noexcept {
        return v;
    }

     [[nodiscard]] ConstantForce mix(const ConstantForce& other) const noexcept {
        return {
            v.x + other.v.x,
            v.y + other.v.y,
            v.z + other.v.z,
            std::max(cutoff_radius, other.cutoff_radius)
        };
    }
};

// Helper to make a dummy particle
static env::impl::Particle make_particle(env::impl::ParticleType type, env::impl::ParticleID id, double mass=1.0, vec3 pos={0,0,0}) {
    return {
    /* id          */ id,
    /* position    */ pos,
    /* velocity    */ vec3{0,0,0},
    /* mass        */ mass,
    /* type        */ type,
    /* state       */ env::impl::Particle::State::ALIVE
    };
}


// Use an environment that supports ConstantForce
using Env = Environment<env::ForcePack<ConstantForce>, env::BoundaryPack<>>;
using IM  = InteractionManager<Env>;
using Info = env::impl::InteractionInfo<IM::force_variant_info_t>; // variant<ConstantForce>


TEST(InteractionManagerTest, EmptyBuild) {
    IM mgr;
    std::vector<Info> info;  // empty

    // empty maps are fine
    EXPECT_NO_THROW(mgr.build(info, {}, {}));
    EXPECT_DOUBLE_EQ(mgr.get_max_cutoff(), 0.0);
}

TEST(InteractionManagerTest, MaxCutoffCalculation) {
    IM mgr;

    // two type-based interactions with cutoffs 1.5 and 2.5
    std::vector<Info> info;
    info.emplace_back(true, std::pair{0, 0}, ConstantForce(1, 1, 1, 1.5));
    info.emplace_back(true, std::pair{1, 1}, ConstantForce(2, 2, 2, 2.5));

    std::unordered_map<ParticleType, env::impl::ParticleType> type_map{{0, 0}, {1, 1}};

    EXPECT_NO_THROW(mgr.build(info, type_map, {}));
    EXPECT_DOUBLE_EQ(mgr.get_max_cutoff(), 2.5);
}

TEST(InteractionManagerTest, TypeBasedLookup) {
    IM mgr;

    std::vector<Info> info;
    info.emplace_back(true, std::pair{0, 0}, ConstantForce(4, 5, 6, -1));
    info.emplace_back(true, std::pair{1, 1}, ConstantForce(1, 2, 3, -1));
    info.emplace_back(true, std::pair{0, 1}, ConstantForce(7, 8, 9, -1));

    std::unordered_map<ParticleType, env::impl::ParticleType> type_map{{0, 0}, {1, 1}};
    mgr.build(info, type_map, {});

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
    IM mgr;

    std::vector<Info> info;
    // Provide a zero type-force for (0,0) so evaluate never hits NullForce
    info.emplace_back(true,  std::pair{0, 0}, ConstantForce(0, 0, 0));
    // one id-based entry for (42,99)
    info.emplace_back(false, std::pair{42, 99}, ConstantForce(7, 8, 9));

    std::unordered_map<ParticleType, env::impl::ParticleType> type_map{{0, 0}};
    std::unordered_map<ParticleID, env::impl::ParticleID> id_map{{42, 0}, {99, 1}};
    mgr.build(info, type_map, id_map);

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
    IM mgr;

    std::vector<Info> info;
    info.emplace_back(true, std::pair{0, 0}, ConstantForce(4, 5, 6, -1));
    info.emplace_back(true, std::pair{1, 1}, ConstantForce(1, 2, 3, -1));

    std::unordered_map<ParticleType, env::impl::ParticleType> type_map{{0, 0}, {1, 1}};
    mgr.build(info, type_map, {});

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
