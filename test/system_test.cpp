#include <gtest/gtest.h>

#include "april/april.hpp"
#include "utils.h"
using namespace april;



TEST(EnvTest, empty_env) {
    Environment e (forces<NoForce>);
    e.set_extent(1,1,1);

    auto sys = build_system(e, DirectSumAoS());

    const auto p = export_particles(sys);
    EXPECT_EQ(p.size(), 0);
}


TEST(EnvTest, one_particle_test) {
    Environment e (forces<LennardJones>);
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10 ));

    e.add_force(LennardJones(3, 5), to_type(0));
    e.set_extent(1,1,1);

    auto sys = build_system(e, DirectSumAoS());
    const auto particles = export_particles(sys);

    EXPECT_EQ(particles.size(), 1);

    const auto p = particles[0];
    EXPECT_TRUE(p.type == 0);
    EXPECT_TRUE(p.id == 0);
    EXPECT_TRUE(p.mass == 10);
    EXPECT_TRUE(p.state == ParticleState::ALIVE);
    EXPECT_TRUE(p.velocity == vec3(1,2,3));
    EXPECT_TRUE(p.position == vec3(3,4,5));
}


TEST(EnvTest, negative_mass_throws) {
    Environment e (forces<NoForce>);
    e.add_particle(make_particle(0, {}, {}, -5));

    e.add_force(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    EXPECT_THROW(build_system(e, DirectSumAoS()), std::invalid_argument);
}



TEST(EnvTest, type_force_missing) {
    Environment e (forces<Gravity>);

    e.add_particle(make_particle(0, {1,2,3}, {0,1,2}, 1, ParticleState::DEAD, -1));
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10, ParticleState::ALIVE, 0));

    e.add_force(Gravity(), between_ids(-1, 0));

    EXPECT_THROW(build_system(e, DirectSumAoS()), std::invalid_argument);
}



TEST(EnvTest, two_particle_force_test) {
    Environment e (forces<Gravity>);

    e.add_particle(make_particle(0, {1,2,3}, {0,1,2}, 1, ParticleState::DEAD, 1));
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10, ParticleState::ALIVE, 0));

    e.add_force(Gravity(), between_ids(1, 0));
    e.add_force(Gravity(), to_type(0));

    auto sys = build_system(e, DirectSumAoS());

    const auto particles = export_particles(sys);
    EXPECT_EQ(particles.size(), 2);

    const auto p1 = particles[0].id == 0? particles[0] : particles[1];
    const auto p2 = particles[0].id == 0? particles[1] : particles[0];

    EXPECT_TRUE(p1.type == 0);
    EXPECT_TRUE(p1.id == 0);
    EXPECT_TRUE(p2.type == 0);
    EXPECT_TRUE(p2.id == 1);
}

TEST(EnvTest, ExtentTooSmallThrows) {
    Environment e (forces<NoForce>);
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE));
    e.add_particle(make_particle(0, {2,0,0}, {}, 1, ParticleState::ALIVE));

    // Set extent too small to cover span=2
    e.set_origin({0,0,0});
    e.set_extent({1,1,1});
    e.add_force(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, OriginOutsideThrows) {
    Environment e (forces<NoForce>);
    // Particles inside [0,1] in each dim
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE));
    e.add_particle(make_particle(0, {1,1,1}, {}, 1, ParticleState::ALIVE));

    // Set origin outside that box
    e.set_origin({2,2,2});
    e.set_extent({2,2,2});
    e.add_force(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, OnlyExtentCentersOrigin) {
    Environment e (forces<NoForce>);
    // Single particle at (3,4,5)
    e.add_particle(make_particle(0, {3,4,5}, {}, 1, ParticleState::ALIVE));

    // Only extent given
    e.set_extent({4,4,4});
    e.add_force(NoForce(), to_type(0));
    const auto sys = build_system(e, DirectSumAoS());

    // bbox_min = (3,4,5), bbox_center = same
    // origin = center - extent/2 = (3,4,5) - (2,2,2) = (1,2,3)
    EXPECT_EQ(sys.domain().origin, vec3(1,2,3));
    EXPECT_EQ(sys.domain().extent, vec3(4,4,4));
}

TEST(EnvTest, OnlyOriginSymmetricExtent) {
    Environment e (forces<NoForce>);
    // Single particle at (3,4,5)
    e.add_particle(make_particle(0, {3,4,5}, {}, 1, ParticleState::ALIVE));

    // Only origin given
    e.set_origin({0,0,0});
    e.add_force(NoForce(), to_type(0));
    e.auto_domain(1);

    const auto sys = build_system(e, DirectSumAoS());


    EXPECT_EQ(sys.domain().origin, vec3(0,0,0));
    EXPECT_EQ(sys.domain().extent, vec3(4,5,6));
}

TEST(EnvTest, AutoOriginExtentDoublesBBox) {
    Environment e (forces<NoForce>);
    // Two particles at (1,2,3) and (3,4,5)
    e.add_particle(make_particle(0, {1,2,3}, {}, 1, ParticleState::ALIVE));
    e.add_particle(make_particle(0, {3,4,5}, {}, 1, ParticleState::ALIVE));

    e.add_force(NoForce(), to_type(0));
    e.auto_domain_factor(1);

    // neither origin nor extent set
    const auto sys = build_system(e, DirectSumAoS());

    // bbox_min = (1,2,3), bbox_max = (3,4,5), bbox_center = (2,3,4), bbox_extent = (2,2,2)
    // extent = bbox_extent * 2 = (4,4,4)
    // origin = center - extent/2 = (2,3,4) - (2,2,2) = (0,1,2)
    EXPECT_EQ(sys.domain().origin, vec3(0,1,2));
    EXPECT_EQ(sys.domain().extent, vec3(4,4,4));
}
