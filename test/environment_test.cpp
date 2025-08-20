#include <gtest/gtest.h>
#include <april/env/environment.h>
#include "april/common.h"

using namespace april;
using namespace april::env;


TEST(EnvTest, empty_env) {
    Environment e;

    e.build();

    auto & p = e.export_particles();
    EXPECT_EQ(p.size(), 0);
}


TEST(EnvTest, one_particle_test) {
    Environment e;

    e.add_particle(Particle{
        .id = PARTICLE_ID_DONT_CARE,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });
    e.add_force_to_type(LennardJones(3,5), 0);
    e.build();

    auto & particles = e.export_particles();
    EXPECT_EQ(particles.size(), 1);

    const impl::Particle & p = particles[0];
    EXPECT_TRUE(p.type == 0);
    EXPECT_TRUE(p.id == 0);
    EXPECT_TRUE(p.mass == 10);
    EXPECT_TRUE(p.state == ParticleState::ALIVE);
    EXPECT_TRUE(p.velocity == vec3(1,2,3));
    EXPECT_TRUE(p.position == vec3(3,4,5));
}


TEST(EnvTest, type_force_missing) {
    Environment e;

    e.add_particle(Particle{
        .id = -1,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add_particle(Particle{
        .id = PARTICLE_ID_DONT_CARE,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add_force_between_ids(InverseSquare(), -1, 0);

    EXPECT_THROW(e.build(), std::invalid_argument);
}



TEST(EnvTest, two_particle_force_test) {
    Environment e;

    e.add_particle(Particle{
        .id = -1,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add_particle(Particle{
        .id = 0,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add_force_between_ids(InverseSquare(), -1, 0);
    e.add_force_to_type(InverseSquare(), 0);

    e.build();

    auto & particles = e.export_particles();
    EXPECT_EQ(particles.size(), 2);

    const impl::Particle & p1 = particles[0].id == 0? particles[0] : particles[1];
    const impl::Particle & p2 = particles[0].id == 0? particles[1] : particles[0];

    EXPECT_TRUE(p1.type == 0);
    EXPECT_TRUE(p1.id == 0);
    EXPECT_TRUE(p2.type == 0);
    EXPECT_TRUE(p2.id == 1);
}


TEST(EnvTest, particle_iterator_test) {
    Environment e;

    e.add_particle(Particle{
        .id = 0,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add_particle(Particle{
        .id = 1,
        .type = 0,
        .position = {3,4,5},
        .velocity = {1,2,3},
        .mass = 10,
        .state = ParticleState::ALIVE,
    });

    e.add_particle(Particle{
        .id = 2,
        .type = 0,
        .position = {1,2,3},
        .velocity = {0,1,2},
        .mass = 1,
        .state = ParticleState::DEAD,
    });

    e.add_force_to_type(NoForce(), 0);

    e.build();

    int i = 0;
    for (const auto & p : e.particles()) {
        EXPECT_TRUE(p.type == 0);
        i++;
    }
    EXPECT_EQ(i, 3);

    i = 0;
    for (const auto & p : e.particles(ParticleState::DEAD)) {
        EXPECT_TRUE(p.mass == 1);
        EXPECT_TRUE(p.state == ParticleState::DEAD);
        i++;
    }
    EXPECT_EQ(i, 2);

    i = 0;
    for (const auto & p : e.particles(ParticleState::ALIVE)) {
        EXPECT_TRUE(p.mass == 10);
        EXPECT_TRUE(p.state == ParticleState::ALIVE);
        i++;
    }
    EXPECT_EQ(i, 1);


    i = 0;
    for (const auto & p : e.particles()) {
        EXPECT_TRUE(p.type == 0);
        i++;
    }
    EXPECT_EQ(i, 3);
}



TEST(EnvTest, ExtentTooSmallThrows) {
    Environment e;
    // Two particles 0 and 2 apart in x
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {0,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {2,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Set extent too small to cover span=2
    e.set_origin({0,0,0});
    e.set_extent({1,1,1});
    e.add_force_to_type(NoForce(), 0);
    EXPECT_THROW(e.build(), std::invalid_argument);
}

TEST(EnvTest, OriginOutsideThrows) {
    Environment e;
    // Particles inside [0,1] in each dim
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {0,0,0}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {1,1,1}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Set origin outside that box
    e.set_origin({2,2,2});
    e.set_extent({2,2,2});
    e.add_force_to_type(NoForce(), 0);
    EXPECT_THROW(e.build(), std::invalid_argument);
}

TEST(EnvTest, OnlyExtentCentersOrigin) {
    Environment e;
    // Single particle at (3,4,5)
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Only extent given
    e.set_extent({4,4,4});
    e.add_force_to_type(NoForce(), 0);
    EXPECT_NO_THROW(e.build());
    const vec3 origin = e.get_origin();
    const vec3 extent = e.get_extent();
    // bbox_min = (3,4,5), bbox_center = same
    // origin = center - extent/2 = (3,4,5) - (2,2,2) = (1,2,3)
    EXPECT_EQ(origin, vec3(1,2,3));
    EXPECT_EQ(extent, vec3(4,4,4));
}

TEST(EnvTest, OnlyOriginSymmetricExtent) {
    Environment e;
    // Single particle at (3,4,5)
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    // Only origin given
    e.set_origin({0,0,0});
    e.add_force_to_type(NoForce(), 0);
    EXPECT_NO_THROW(e.build());
    const vec3 origin = e.get_origin();
    const vec3 extent = e.get_extent();
    // bbox_center = (3,4,5), opposite = origin + 2*(center-origin) = 2*center = (6,8,10)
    // extent = abs(opposite-origin) = (6,8,10)
    EXPECT_EQ(origin, vec3(0,0,0));
    EXPECT_EQ(extent, vec3(6,8,10));
}

TEST(EnvTest, AutoOriginExtentDoublesBBox) {
    Environment e;
    // Two particles at (1,2,3) and (3,4,5)
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {1,2,3}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add_particle({.id = PARTICLE_ID_DONT_CARE, .type = 0, .position = {3,4,5}, .velocity = {0,0,0}, .mass = 1, .state = ParticleState::ALIVE});
    e.add_force_to_type(NoForce(), 0);
    // neither origin nor extent set
    EXPECT_NO_THROW(e.build());
    const vec3 origin = e.get_origin();
    const vec3 extent = e.get_extent();
    // bbox_min = (1,2,3), bbox_max = (3,4,5), bbox_center = (2,3,4), bbox_extent = (2,2,2)
    // extent = bbox_extent * 2 = (4,4,4)
    // origin = center - extent/2 = (2,3,4) - (2,2,2) = (0,1,2)
    EXPECT_EQ(origin, vec3(0,1,2));
    EXPECT_EQ(extent, vec3(4,4,4));
}
