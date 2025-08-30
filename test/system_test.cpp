#include <gtest/gtest.h>
#include <april/env/environment.h>
#include "april/common.h"
#include <april/core/system.h>
#include <april/containers/direct_sum.h>

using namespace april;
using namespace april::env;
using namespace april::core;
using namespace april::cont;


TEST(EnvTest, empty_env) {
    Environment e (forces<NoForce>);

    auto sys = build_system(e, DirectSum());

    const auto p = sys.export_particles();
    EXPECT_EQ(p.size(), 0);
}


TEST(EnvTest, one_particle_test) {
    Environment e (forces<LennardJones>);
    e.add(Particle{
        .id = PARTICLE_ID_DONT_CARE,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });
    e.add_force(LennardJones(3, 5), to_type(0));

    auto sys = build_system(e, DirectSum());
    auto particles = sys.export_particles();

    EXPECT_EQ(particles.size(), 1);

    const env::impl::ParticleView p = particles[0];
    EXPECT_TRUE(p.type == 0);
    EXPECT_TRUE(p.id == 0);
    EXPECT_TRUE(p.mass == 10);
    EXPECT_TRUE(p.state == ParticleState::ALIVE);
    EXPECT_TRUE(p.velocity == vec3(1,2,3));
    EXPECT_TRUE(p.position == vec3(3,4,5));
}


TEST(EnvTest, type_force_missing) {
    Environment e (forces<InverseSquare>);

    e.add(Particle{
        .id = -1,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add(Particle{
        .id = PARTICLE_ID_DONT_CARE,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add_force(InverseSquare(), between_ids(-1, 0));

    EXPECT_THROW(build_system(e, DirectSum()), std::invalid_argument);
}



TEST(EnvTest, two_particle_force_test) {
    Environment e (forces<InverseSquare>);

    e.add(Particle{
        .id = -1,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add(Particle{
        .id = 0,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add_force(InverseSquare(), between_ids(-1, 0));
    e.add_force(InverseSquare(), to_type(0));


    auto sys = build_system(e, DirectSum());

    auto particles = sys.export_particles();
    EXPECT_EQ(particles.size(), 2);

    const env::impl::ParticleView & p1 = particles[0].id == 0? particles[0] : particles[1];
    const env::impl::ParticleView & p2 = particles[0].id == 0? particles[1] : particles[0];

    EXPECT_TRUE(p1.type == 0);
    EXPECT_TRUE(p1.id == 0);
    EXPECT_TRUE(p2.type == 0);
    EXPECT_TRUE(p2.id == 1);
}


TEST(EnvTest, particle_iterator_test) {
    Environment e (forces<NoForce>);

    e.add(Particle{
        .id = 0,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add(Particle{
        .id = 1,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add(Particle{
        .id = 2,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add_force(NoForce(), to_type(0));

    auto sys = build_system(e, DirectSum());

    int i = 0;
    for (const auto & p : sys.export_particles()) {
        EXPECT_TRUE(p.type == 0);
        i++;
    }
    EXPECT_EQ(i, 3);

    i = 0;
    for (const auto & p : sys.export_particles(ParticleState::DEAD)) {
        EXPECT_TRUE(p.mass == 1);
        EXPECT_TRUE(p.state == ParticleState::DEAD);
        i++;
    }
    EXPECT_EQ(i, 2);

    i = 0;
    for (const auto & p : sys.export_particles(ParticleState::ALIVE)) {
        EXPECT_TRUE(p.mass == 10);
        EXPECT_TRUE(p.state == ParticleState::ALIVE);
        i++;
    }
    EXPECT_EQ(i, 1);


    i = 0;
    for (const auto & p : sys.export_particles()) {
        EXPECT_TRUE(p.type == 0);
        i++;
    }
    EXPECT_EQ(i, 3);
}



TEST(EnvTest, ExtentTooSmallThrows) {
    Environment e (forces<NoForce>);
    // Two particles 0 and 2 apart in x
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {0,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {2,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Set extent too small to cover span=2
    e.set_origin({0,0,0});
    e.set_extent({1,1,1});
    e.add_force(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, DirectSum()), std::invalid_argument);
}

TEST(EnvTest, OriginOutsideThrows) {
    Environment e (forces<NoForce>);
    // Particles inside [0,1] in each dim
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {0,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {1,1,1}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Set origin outside that box
    e.set_origin({2,2,2});
    e.set_extent({2,2,2});
    e.add_force(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, DirectSum()), std::invalid_argument);
}

TEST(EnvTest, OnlyExtentCentersOrigin) {
    Environment e (forces<NoForce>);
    // Single particle at (3,4,5)
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Only extent given
    e.set_extent({4,4,4});
    e.add_force(NoForce(), to_type(0));
    auto sys = build_system(e, DirectSum());
    const vec3 origin = sys.domain.origin;
    const vec3 extent = sys.domain.extent;
    // bbox_min = (3,4,5), bbox_center = same
    // origin = center - extent/2 = (3,4,5) - (2,2,2) = (1,2,3)
    EXPECT_EQ(origin, vec3(1,2,3));
    EXPECT_EQ(extent, vec3(4,4,4));
}

TEST(EnvTest, OnlyOriginSymmetricExtent) {
    Environment e (forces<NoForce>);
    // Single particle at (3,4,5)
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Only origin given
    e.set_origin({0,0,0});
    e.add_force(NoForce(), to_type(0));
    auto sys = build_system(e, DirectSum());

    const vec3 origin = sys.domain.origin;
    const vec3 extent = sys.domain.extent;
    // bbox_center = (3,4,5), opposite = origin + 2*(center-origin) = 2*center = (6,8,10)
    // extent = abs(opposite-origin) = (6,8,10)
    EXPECT_EQ(origin, vec3(0,0,0));
    EXPECT_EQ(extent, vec3(6,8,10));
}

TEST(EnvTest, AutoOriginExtentDoublesBBox) {
    Environment e (forces<NoForce>);
    // Two particles at (1,2,3) and (3,4,5)
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {1,2,3}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add_force(NoForce(), to_type(0));
    // neither origin nor extent set
    auto sys = build_system(e, DirectSum());

    const vec3 origin = sys.domain.origin;
    const vec3 extent = sys.domain.extent;
    // bbox_min = (1,2,3), bbox_max = (3,4,5), bbox_center = (2,3,4), bbox_extent = (2,2,2)
    // extent = bbox_extent * 2 = (4,4,4)
    // origin = center - extent/2 = (2,3,4) - (2,2,2) = (0,1,2)
    EXPECT_EQ(origin, vec3(0,1,2));
    EXPECT_EQ(extent, vec3(4,4,4));
}
