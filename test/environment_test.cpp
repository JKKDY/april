#include <gtest/gtest.h>
#include <april/env/environment.h>

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
        .id = PARTICLE_ID_UNDEFINED,
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
    EXPECT_TRUE(p.velocity == april::utils::Vec3(1,2,3));
    EXPECT_TRUE(p.position == april::utils::Vec3(3,4,5));
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
        .id = PARTICLE_ID_UNDEFINED,
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
};