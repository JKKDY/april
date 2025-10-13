#include <gtest/gtest.h>
#include "april/april.h"

using namespace april;


struct TouchSpy final : Boundary {
	using PID = env::internal::ParticleID;

	// thickness >= 0 → inside slab; < 0 → outside half-space
	explicit TouchSpy(double thickness, std::vector<PID>* sink)
	: Boundary(thickness, false, false, false), sink(sink) {}

	void apply(env::internal::Particle& p, const env::Box &, boundary::Face) const noexcept {
		if (sink) sink->push_back(p.id);
	}

private:
	std::vector<PID>* sink;
};


template <class ContainerT>
class BoundaryTestT : public testing::Test {};

using ContainerTypes = testing::Types<DirectSum, LinkedCells>;
TYPED_TEST_SUITE(BoundaryTestT, ContainerTypes);

TYPED_TEST(BoundaryTestT, InsideSlab_XMinus_AppliesOnlyToSlabParticles) {
	// Domain: [0,10]^3
	Environment env (forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});

	// Particles: one in the X- slab [0,1], one outside it
	env.add({.id=0, .type=0, .position={0.4,5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE});
	env.add({.id=1, .type=0, .position={2.0,5,5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE});

	env.add_force(NoForce{}, to_type(0));

	// External sinks for each face
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// Set per-face spies (thickness 1.0 on X-, 0 elsewhere)
	env.set_boundaries(std::array{
		TouchSpy{1.0, &sxm}, // X-
		TouchSpy{0.0, &sxp}, // X+
		TouchSpy{0.0, &sym}, // Y-
		TouchSpy{0.0, &syp}, // Y+
		TouchSpy{0.0, &szm}, // Z-
		TouchSpy{0.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;

	auto sys = build_system(env, TypeParam(), &mappings);
	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id0 = mappings.usr_ids_to_impl_ids.at(0);
	ASSERT_EQ(sxm.size(), 1u);
	EXPECT_EQ(sxm[0], id0);
	EXPECT_TRUE(sxp.empty());
	EXPECT_TRUE(sym.empty());
	EXPECT_TRUE(syp.empty());
	EXPECT_TRUE(szm.empty());
	EXPECT_TRUE(szp.empty());
}


TYPED_TEST(BoundaryTestT, OutsideHalfspace_XPlus_TouchesOnlyActualExiters) {
	Environment env (forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// p0: crosses X+ this step
	env.add({.id=0, .type=0, .position={9.5,5,5}, .velocity={+2,0,0}, .mass=1, .state=ParticleState::ALIVE});
	// p1: exits via Y+ instead
	env.add({.id=1, .type=0, .position={5.0,9.5,5}, .velocity={0,+2,0}, .mass=1, .state=ParticleState::ALIVE});

	// external sinks
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// boundaries: only X+ and Y+ are outside (thickness < 0)
	env.set_boundaries(std::array{
		TouchSpy{ 0.0, &sxm}, // X-
		TouchSpy{-1.0, &sxp}, // X+
		TouchSpy{ 0.0, &sym}, // Y-
		TouchSpy{-1.0, &syp}, // Y+
		TouchSpy{ 0.0, &szm}, // Z-
		TouchSpy{ 0.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// mark old positions behind where they will move from
	for (auto pid = sys.index_start(); pid < sys.index_end(); pid++) {
		env::internal::Particle & p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // simulate one step
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id0 = mappings.usr_ids_to_impl_ids.at(0);
	auto id1 = mappings.usr_ids_to_impl_ids.at(1);

	ASSERT_EQ(sxp.size(), 1u);
	EXPECT_EQ(sxp[0], id0);    // p0 → X+
	ASSERT_EQ(syp.size(), 1u);
	EXPECT_EQ(syp[0], id1);    // p1 → Y+
	EXPECT_TRUE(sxm.empty());
	EXPECT_TRUE(sym.empty());
	EXPECT_TRUE(szm.empty());
	EXPECT_TRUE(szp.empty());
}


TYPED_TEST(BoundaryTestT, CornerExit_TriggersRelevantFaces) {
	Environment env (forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// One particle moving diagonally out through the X+Y+ edge
	env.add({.id=42, .type=0, .position={9.7,9.7,5}, .velocity={+1,+1,0}, .mass=1, .state=ParticleState::ALIVE});

	// sinks
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// outside boundaries on X+ and Y+ only
	env.set_boundaries(std::array{
		TouchSpy{ 0.0, &sxm}, // X-
		TouchSpy{-1.0, &sxp}, // X+
		TouchSpy{ 0.0, &sym}, // Y-
		TouchSpy{-1.0, &syp}, // Y+
		TouchSpy{ 0.0, &szm}, // Z-
		TouchSpy{ 0.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// mark old positions behind where they will move from
	for (auto pid = sys.index_start(); pid < sys.index_end(); ++pid) {
		env::internal::Particle & p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // simulate one step
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id42 = mappings.usr_ids_to_impl_ids.at(42);

	bool x_hit = std::ranges::find(sxp, id42) != sxp.end();
	bool y_hit = std::ranges::find(syp, id42) != syp.end();

	EXPECT_TRUE(x_hit || y_hit)
		<< "Particle should trigger at least one of X+ or Y+ faces at the corner";
	EXPECT_TRUE(sxm.empty());
	EXPECT_TRUE(sym.empty());
	EXPECT_TRUE(szm.empty());
	EXPECT_TRUE(szp.empty());
}


TYPED_TEST(BoundaryTestT, InsideCorner_TouchesAllOverlappingFaces) {
	// Domain [0,10]^3
	Environment env(forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Place a particle inside the corner region (X−, Y−, Z−)
	// so it's within all three inside slabs of thickness=1.
	env.add({.id=0, .type=0, .position={0.5,0.5,0.5},
			 .velocity={}, .mass=1, .state=ParticleState::ALIVE});

	// External sinks for each face
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// Inside slabs on all faces (thickness=1)
	env.set_boundaries(std::array{
		TouchSpy{1.0, &sxm}, // X−
		TouchSpy{1.0, &sxp}, // X+
		TouchSpy{1.0, &sym}, // Y−
		TouchSpy{1.0, &syp}, // Y+
		TouchSpy{1.0, &szm}, // Z−
		TouchSpy{1.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);
	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id = mappings.usr_ids_to_impl_ids.at(0);

	// The particle should be touched by the three minus faces.
	EXPECT_EQ(sxm, std::vector{id});
	EXPECT_EQ(sym, std::vector{id});
	EXPECT_EQ(szm, std::vector{id});

	// No contact with opposite (+) faces.
	EXPECT_TRUE(sxp.empty());
	EXPECT_TRUE(syp.empty());
	EXPECT_TRUE(szp.empty());
}


TYPED_TEST(BoundaryTestT, NearCornerExit_TriggersCorrectFace) {
	Environment env (forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// One particle moving diagonally through the +z face but ending up in the overlap region of the +x,+y,+z boundaries
	env.add({.id=42, .type=0, .position={9.7,9.7,9.8}, .velocity={1,1,1}, .mass=1, .state=ParticleState::ALIVE});

	// sinks
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// outside boundaries on X+ and Y+ only
	env.set_boundaries(std::array{
		TouchSpy{-1.0, &sxm}, // X-
		TouchSpy{-1.0, &sxp}, // X+
		TouchSpy{-1.0, &sym}, // Y-
		TouchSpy{-1.0, &syp}, // Y+
		TouchSpy{-1.0, &szm}, // Z-
		TouchSpy{-1.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// mark old positions behind where they will move from
	for (auto pid = sys.index_start(); pid < sys.index_end(); pid++) {
		env::internal::Particle & p = sys.get_particle_by_index(pid);
		p.old_position = p.position;
		p.position = p.old_position + p.velocity; // simulate one step
	}

	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	auto id = mappings.usr_ids_to_impl_ids.at(42);

	EXPECT_TRUE(sxp.empty());
	EXPECT_TRUE(sxm.empty());
	EXPECT_TRUE(syp.empty());
	EXPECT_TRUE(sym.empty());
	EXPECT_TRUE(szm.empty());
	EXPECT_EQ(szp.size(), 1);
	EXPECT_EQ(szp[0], id);
}


TYPED_TEST(BoundaryTestT, InsideSlab_AllFaces_OneParticleEach) {
	Environment env(forces<NoForce>, boundary::boundaries<TouchSpy>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Six particles, one for each face region.
	// Positions are clearly inside their slabs (thickness = 1)
	env.add({.id=0, .type=0, .position={0.5, 5, 5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // X−
	env.add({.id=1, .type=0, .position={9.5, 5, 5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // X+
	env.add({.id=2, .type=0, .position={5, 0.5, 5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Y−
	env.add({.id=3, .type=0, .position={5, 9.5, 5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Y+
	env.add({.id=4, .type=0, .position={5, 5, 0.5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Z−
	env.add({.id=5, .type=0, .position={5, 5, 9.5}, .velocity={}, .mass=1, .state=ParticleState::ALIVE}); // Z+

	// External sinks for each face
	std::vector<env::internal::ParticleID> sxm, sxp, sym, syp, szm, szp;

	// Assign spies: thickness 1.0 for all faces
	env.set_boundaries(std::array{
		TouchSpy{1.0, &sxm}, // X−
		TouchSpy{1.0, &sxp}, // X+
		TouchSpy{1.0, &sym}, // Y−
		TouchSpy{1.0, &syp}, // Y+
		TouchSpy{1.0, &szm}, // Z−
		TouchSpy{1.0, &szp}  // Z+
	});

	UserToInternalMappings mappings;
	auto sys = build_system(env, TypeParam(), &mappings);
	sys.register_all_particle_movements();
	sys.apply_boundary_conditions();

	// Map user IDs → internal IDs for verification
	auto get_id = [&](int uid) { return mappings.usr_ids_to_impl_ids.at(uid); };

	EXPECT_EQ(sxm, std::vector{get_id(0)});
	EXPECT_EQ(sxp, std::vector{get_id(1)});
	EXPECT_EQ(sym, std::vector{get_id(2)});
	EXPECT_EQ(syp, std::vector{get_id(3)});
	EXPECT_EQ(szm, std::vector{get_id(4)});
	EXPECT_EQ(szp, std::vector{get_id(5)});
}





