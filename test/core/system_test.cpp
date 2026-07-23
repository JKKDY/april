#include <gtest/gtest.h>

#include "april/april.hpp"
#include "utils.h"
using namespace april;



TEST(EnvTest, empty_env) {
    Environment e (forces<NoForce>);
    e.set_extent(1,1,1);

    auto sys = build_system(e, container::DirectSumAoS());

    const auto p = export_particles(sys);
    EXPECT_EQ(p.size(), 0);
}


TEST(EnvTest, one_particle_test) {
    Environment e (forces<LennardJones>);
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10 ));

    e.add_interaction(LennardJones(3, 5), to_type(0));
    e.set_extent(1,1,1);

    auto sys = build_system(e, container::DirectSumAoS());
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

    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}



TEST(EnvTest, type_force_missing) {
    Environment e (forces<Gravity>);

    e.add_particle(make_particle(0, {1,2,3}, {0,1,2}, 1, ParticleState::DEAD, -1));
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10, ParticleState::ALIVE, 0));

    e.add_interaction(Gravity(), between_ids(-1, 0));

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}



TEST(EnvTest, two_particle_force_test) {
    Environment e (forces<Gravity>);

    e.add_particle(make_particle(0, {1,2,3}, {0,1,2}, 1, ParticleState::DEAD, 1));
    e.add_particle(make_particle(0, {3,4,5}, {1,2,3}, 10, ParticleState::ALIVE, 0));

    e.add_interaction(Gravity(), between_ids(1, 0));
    e.add_interaction(Gravity(), to_type(0));

    auto sys = build_system(e, container::DirectSumAoS());

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
    e.add_interaction(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, OriginOutsideThrows) {
    Environment e (forces<NoForce>);
    // Particles inside [0,1] in each dim
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE));
    e.add_particle(make_particle(0, {1,1,1}, {}, 1, ParticleState::ALIVE));

    // Set origin outside that box
    e.set_origin({2,2,2});
    e.set_extent({2,2,2});
    e.add_interaction(NoForce(), to_type(0));
    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, OnlyExtentCentersOrigin) {
    Environment e (forces<NoForce>);
    // Single particle at (3,4,5)
    e.add_particle(make_particle(0, {3,4,5}, {}, 1, ParticleState::ALIVE));

    // Only extent given
    e.set_extent({4,4,4});
    e.add_interaction(NoForce(), to_type(0));
    const auto sys = build_system(e, container::DirectSumAoS());

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
    e.add_interaction(NoForce(), to_type(0));
    e.set_domain_padding(1);

    const auto sys = build_system(e, container::DirectSumAoS());


    EXPECT_EQ(sys.domain().origin, vec3(0,0,0));
    EXPECT_EQ(sys.domain().extent, vec3(4,5,6));
}

TEST(EnvTest, DomainPaddingFactorAddsRelativePaddingPerSide) {
    Environment e(forces<NoForce>);

    e.add_particle(make_particle(0, {1,2,3}, {}, 1, ParticleState::ALIVE));
    e.add_particle(make_particle(0, {3,4,5}, {}, 1, ParticleState::ALIVE));

    e.add_interaction(NoForce(), to_type(0));
    e.set_domain_padding_factor(2);

    // Neither origin nor extent set.
    const auto sys = build_system(e, container::DirectSumAoS());

    // bbox_min    = (1,2,3)
    // bbox_max    = (3,4,5)
    // bbox_center = (2,3,4)
    // bbox_extent = (2,2,2)
    //
    // padding_per_side = bbox_extent * 2 = (4,4,4)
    // extent           = bbox_extent + 2 * padding_per_side = (10,10,10)
    // origin           = center - extent/2 = (2,3,4) - (5,5,5) = (-3,-2,-1)
    EXPECT_EQ(sys.domain().origin, vec3(-3,-2,-1));
    EXPECT_EQ(sys.domain().extent, vec3(10,10,10));
}


TEST(EnvTest, IdentityMappingForDenseInput) {
    // Verify that dense, sorted, non-interacting input
    // results in an identity mapping (User ID == System ID).
    Environment e(forces<NoForce>);

    // Input: IDs 0, 1, 2 | Types 0, 1
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
    e.add_particle(make_particle(1, {1,1,1}, {}, 1, ParticleState::ALIVE, 1));
    e.add_particle(make_particle(0, {2,2,2}, {}, 1, ParticleState::ALIVE, 2));

    // Self-interactions only (no ID-to-ID to trigger reordering)
    e.add_interaction(NoForce(), to_type(0));
    e.add_interaction(NoForce(), to_type(1));
    e.set_extent(10, 10, 10);

    BuildInfo info;
    auto sys = build_system(e, container::DirectSumAoS(), ExecutionConfig(), &info);

    // Verify Type Identity
    EXPECT_EQ(info.type_map.at(0), 0);
    EXPECT_EQ(info.type_map.at(1), 1);

    // Verify ID Identity
    EXPECT_EQ(info.id_map.at(0), 0);
    EXPECT_EQ(info.id_map.at(1), 1);
    EXPECT_EQ(info.id_map.at(2), 2);
}


TEST(EnvTest, StableMappingAndInteractionPrioritization) {
    // We want to verify:
    // 1. Interacting IDs are at the front of the system memory.
    // 2. Relative order of IDs is preserved (Stable mapping).
    // 3. User types map to dense indices in ascending order.

    Environment e(forces<NoForce>);

    // Create sparse, out-of-order IDs and types
    // Non-interacting IDs: 100, 200
    // Interacting IDs: 10, 50
    // Types: 5, 2
    e.add_particle(make_particle(5, {0,0,0}, {}, 1, ParticleState::ALIVE, 200));
    e.add_particle(make_particle(2, {1,1,1}, {}, 1, ParticleState::ALIVE, 50));
    e.add_particle(make_particle(5, {2,2,2}, {}, 1, ParticleState::ALIVE, 10));
    e.add_particle(make_particle(2, {3,3,3}, {}, 1, ParticleState::ALIVE, 100));

    // Define interactions to trigger prioritization
    e.add_interaction(NoForce(), between_ids(10, 50));
    // Self-interactions required by validation
    e.add_interaction(NoForce(), to_type(2));
    e.add_interaction(NoForce(), to_type(5));

    e.set_extent(10, 10, 10);

    BuildInfo info;
    auto sys = build_system(e, container::DirectSumAoS(), ExecutionConfig(), &info);

    // Check Type Mapping (Ascending User Type -> Dense System Index)
    // User Type 2 -> System Type 0
    // User Type 5 -> System Type 1
    EXPECT_EQ(info.type_map.at(2), 0);
    EXPECT_EQ(info.type_map.at(5), 1);

    // Check ID Mapping (Interacting IDs at front, then others, all stable)
    // Expected order in memory: 10, 50, 100, 200
    // (10, 50 are interacting, 100, 200 are not. Within groups, they stay ascending)
    EXPECT_EQ(info.id_map.at(10),  0);
    EXPECT_EQ(info.id_map.at(50),  1);
    EXPECT_EQ(info.id_map.at(100), 2);
    EXPECT_EQ(info.id_map.at(200), 3);
}


TEST(EnvTest, SparseAndAutoIDMix) {
    Environment e(forces<NoForce>);
    // User IDs: 1 and 10. Missing IDs for others.
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 1));
    e.add_particle(make_particle(0, {1,1,1}, {}, 1, ParticleState::ALIVE, std::nullopt));
    e.add_particle(make_particle(0, {2,2,2}, {}, 1, ParticleState::ALIVE, 10));
    e.add_particle(make_particle(0, {3,3,3}, {}, 1, ParticleState::ALIVE, std::nullopt));

    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(10,10,10);

    BuildInfo info;
    build_system(e, container::DirectSumAoS(), ExecutionConfig(), &info);

    // Should find 4 unique mappings. Auto-ids should skip 1 and 10.
    EXPECT_EQ(info.id_map.size(), 4);
    EXPECT_TRUE(info.id_map.contains(1));
    EXPECT_TRUE(info.id_map.contains(10));
    // Auto-assigned should be 0 and 2 (first available non-colliding IDs)
    EXPECT_TRUE(info.id_map.contains(0));
    EXPECT_TRUE(info.id_map.contains(2));
}


TEST(EnvTest, MissingSelfInteractionThrows) {
    Environment e(forces<NoForce>);
    e.add_particle(make_particle(0, {0,0,0}, {}, 1));
    e.add_particle(make_particle(1, {1,1,1}, {}, 1));

    // Interaction between 0 and 1 exists, but 0-0 and 1-1 are missing
    e.add_interaction(NoForce(), between_types(0, 1));

    e.set_extent(10, 10, 10);
    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, InvalidStateThrows) {
    Environment e(forces<NoForce>);
    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(10,10,10);

    // Case A: INVALID sentinel
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::INVALID, 0));
    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);

    // Case B: Undefined bits
    Environment e2(forces<NoForce>);
    e2.add_interaction(NoForce(), to_type(0));
    e2.set_extent(10,10,10);
    e2.add_particle(make_particle(0, {0,0,0}, {}, 1, static_cast<ParticleState>(0b10101010), 1));
    EXPECT_THROW(build_system(e2, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, IDInteractionPriority) {
    Environment e(forces<NoForce>);
    // IDs 0, 1, 2. All are Type 0.
    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
    e.add_particle(make_particle(0, {1,1,1}, {}, 1, ParticleState::ALIVE, 1));
    e.add_particle(make_particle(0, {2,2,2}, {}, 1, ParticleState::ALIVE, 2));

    e.add_interaction(NoForce(), to_type(0));
    // Force between 1 and 2 should move them to system indices 0 and 1
    e.add_interaction(NoForce(), between_ids(1, 2));

    e.set_extent(10,10,10);

    BuildInfo info;
    build_system(e, container::DirectSumAoS(), ExecutionConfig(), &info);

    // User ID 0 was not in an ID-interaction, should be pushed to the back (index 2)
    EXPECT_EQ(info.id_map.at(1), 0);
    EXPECT_EQ(info.id_map.at(2), 1);
    EXPECT_EQ(info.id_map.at(0), 2);
}



TEST(EnvTest, AbsoluteMarginExpansion) {
    Environment e(forces<NoForce>);
    e.add_particle(vec3(0,0,0), vec3(0,0,0), 1.0);
    e.add_interaction(NoForce(), to_type(0));

    // Uses your API: auto_domain(double) sets absolute margin
    e.set_domain_padding(1.0);

    const auto sys = build_system(e, container::DirectSumAoS());

    // Particle at [0,0,0], abs margin 1.0 -> Domain [-1, 1]
    EXPECT_EQ(sys.domain().origin, vec3(-1,-1,-1));
    EXPECT_EQ(sys.domain().extent, vec3(2,2,2));
}

TEST(EnvTest, ZeroExtentThrows) {
    Environment e(forces<NoForce>);
    e.add_particle(vec3(0,0,0), vec3(0,0,0), 1.0);
    e.add_interaction(NoForce(), to_type(0));

    // Use chaining to set both to zero
    e.with_domain_padding(0.0).set_domain_padding_factor(vec3(0,0,0));

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::logic_error);
}

TEST(EnvTest, NegativeMarginsThrow) {
    Environment e(forces<NoForce>);
    e.add_particle(vec3(0,0,0), vec3(0,0,0), 1.0);
    e.add_interaction(NoForce(), to_type(0));

    // Passing negative to the vec3 overload
    e.set_domain_padding(vec3(-1, 0, 0));
    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::logic_error);
}

TEST(EnvTest, MarginPriorityMax) {
    Environment e(forces<NoForce>);
    // BBox extent is 10.0
    e.add_particle(vec3(0,0,0), vec3(0,0,0), 1.0);
    e.add_particle(vec3(10,0,0), vec3(0,0,0), 1.0);
    e.add_interaction(NoForce(), to_type(0));

    // Factor 0.1 of 10 = 1.0
    // Absolute margin = 5.0
    e.set_domain_padding_factor(vec3(0.1, 0.1, 0.1));
    e.set_domain_padding(5.0);

    const auto sys = build_system(e, container::DirectSumAoS());

    // Absolute margin (5) > Factor margin (1)
    // origin.x = bbox_min(0) - 5 = -5
    EXPECT_EQ(sys.box().min.x, -5.0);
}


TEST(EnvTest, DuplicateIDThrows) {
    Environment e(forces<NoForce>);
    Particle p = make_particle(0, {0,0,0}, {0,0,0}, 1.0, ParticleState::ALIVE, 10);

    e.add_particle(p);
    // Adding the same ID again must throw immediately
    EXPECT_THROW(e.add_particle(p), std::invalid_argument);
}


TEST(EnvTest, NonExistingTypeInteractionThrows) {
    Environment e(forces<NoForce>);
    e.add_particle(vec3(0), vec3(0), 1.0, 0); // Only Type 0 exists

    // Force assigned to Type 999
    e.add_interaction(NoForce(), to_type(999));
    e.set_extent(1,1,1);

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}


TEST(EnvTest, NonExistingIDInteractionThrows) {
    Environment e(forces<NoForce>);
    e.add_particle(make_particle(0, {0,0,0}, {0,0,0}, 1.0, ParticleState::ALIVE, 1));

    // Force between existing ID 1 and non-existing ID 999
    e.add_interaction(NoForce(), between_ids(1, 999));
    e.set_extent(1,1,1);

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, ZeroStateThrows) {
    Environment e(forces<NoForce>);
    // Manual cast to bypass ALIVE default
    e.add_particle(make_particle(0, {0,0,0}, {0,0,0}, 1.0, static_cast<ParticleState>(0), 1));

    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}


TEST(EnvTest, DefaultToOpenBoundary) {
    Environment e(forces<NoForce>);
    e.add_particle(vec3(0), vec3(0), 1.0);
    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    // No boundaries set manually
    EXPECT_NO_THROW(build_system(e, container::DirectSumAoS()));
}

TEST(EnvTest, AsymmetricPeriodicThrows) {
    // Assuming PeriodicBoundary is available in your traits
    Environment e(forces<NoForce>, boundaries<PeriodicBoundary, OpenBoundary>);
    e.add_particle(vec3(0.5), vec3(0), 1.0);
    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    // Set only X- to Periodic; X+ defaults to Open
    e.set_boundary(PeriodicBoundary(), DomainFace::XMinus);

    EXPECT_THROW(build_system(e, container::DirectSumAoS()), std::invalid_argument);
}

TEST(EnvTest, PeriodicityFlagsPropagate) {
    Environment e(forces<NoForce>, boundaries<PeriodicBoundary>);
    e.add_particle(vec3(0.5), vec3(0), 1.0);
    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    // Symmetrically set X axis to Periodic
    e.set_boundary(PeriodicBoundary(), DomainFace::XMinus);
    e.set_boundary(PeriodicBoundary(), DomainFace::XPlus);

    auto sys = build_system(e, container::DirectSumAoS());
}

TEST(EnvTest, MixedBoundaries) {
    Environment e(forces<NoForce>, boundaries<PeriodicBoundary, ReflectiveBoundary>);
    e.add_particle(vec3(0.5), vec3(0), 1.0);
    e.add_interaction(NoForce(), to_type(0));
    e.set_extent(1,1,1);

    // x-axis: Periodic
    e.set_boundary(PeriodicBoundary(), DomainFace::XMinus);
    e.set_boundary(PeriodicBoundary(), DomainFace::XPlus);

    // y-axis: Reflective
    e.set_boundary(ReflectiveBoundary(), DomainFace::YMinus);
    e.set_boundary(ReflectiveBoundary(), DomainFace::YPlus);

    EXPECT_NO_THROW(build_system(e, container::DirectSumAoS()));
}

