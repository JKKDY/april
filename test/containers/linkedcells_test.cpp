#include <gtest/gtest.h>
#include <gmock/gmock.h>


using testing::AnyOf;
using testing::Eq;


#include "april/containers/linked_cells.hpp"


#include "orbit_monitor.h"
#include "constant_force.h"
#include "utils.h"

using namespace april;





// Execution Strategy Wrapper
template<ParallelPolicy P, VectorPolicy V = VectorPolicy::Auto>
struct CustomExecConfig : RunTimeConfig<>, CompileTimeConfig<P, V> {};

// Ordering Policies
struct OrderDefault {
    static auto&& apply(auto&& c) { return std::forward<decltype(c)>(c); }
};

struct OrderHilbert {
    static auto&& apply(auto&& c) { return c.with_cell_ordering(hilbert_order); }
};

template <typename ContainerT, typename OrderingT, ParallelPolicy P, VectorPolicy V>
struct TestConfig {
    using Container = ContainerT;
    using ExecConfig = CustomExecConfig<P, V>;

    // Returns the configured container by value (ready for build_system)
    static auto create_container(double cell_size) {
        auto c = ContainerT{};

        c.with_abs_cell_size(cell_size)
         .with_skin_factor(0.0)
         .with_block_size(2);

        return OrderingT::apply(std::move(c));
    }

	static auto create_container() {
    	auto c = ContainerT{};
    	c.with_skin_factor(0.0);
    	return OrderingT::apply(std::move(c));
    }

    // Returns a valid execution config
    static auto create_exec() {
        return ExecConfig{};
    }
};

// 4. Optimized Test Matrix (The "Five-Point" Coverage)
using Matrix = testing::Types<
    TestConfig<LinkedCells<Layout::AoS>, OrderDefault, ParallelPolicy::Serial, VectorPolicy::Scalar>,

    TestConfig<LinkedCells<Layout::SoA>, OrderDefault, ParallelPolicy::Serial, VectorPolicy::Auto>,

    TestConfig<LinkedCells<Layout::SoA>, OrderDefault, ParallelPolicy::Threaded, VectorPolicy::Auto>,

    TestConfig<LinkedCells<Layout::AoSoA<8>>, OrderDefault, ParallelPolicy::Serial, VectorPolicy::Scalar>,

	TestConfig<LinkedCells<Layout::AoSoA<8>>, OrderDefault, ParallelPolicy::Serial, VectorPolicy::Auto>,

    TestConfig<LinkedCells<Layout::AoSoA<32>>, OrderHilbert, ParallelPolicy::Threaded, VectorPolicy::Auto>
>;

template <typename T>
class LinkedCellsTest : public testing::Test {};

TYPED_TEST_SUITE(LinkedCellsTest, Matrix);



TYPED_TEST(LinkedCellsTest, SingleParticle_NoForce) {
    Environment e (forces<NoForce>);
	e.add_particle(make_particle(0, {1,2,3}, {}, 1, ParticleState::ALIVE, 0));
	e.add_force(NoForce(), to_type(0));
	e.set_extent({4,4,4});

	auto sys = build_system(e, TypeParam::create_container(4), TypeParam::create_exec());
    sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 1u);
    EXPECT_EQ(out[0].force, vec3(0,0,0));
}

TYPED_TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_SameCell) {
    Environment e(forces<ConstantForce>);
	e.set_extent({2,2,2});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1.5,0,0}, {}, 2, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, TypeParam::create_container(2), TypeParam::create_exec());
	sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, p2.force);
}

TYPED_TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,1});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1.5,0,0}, {}, 2, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, TypeParam::create_container(1), TypeParam::create_exec());
	sys.update_forces();

    auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	EXPECT_THAT(p1.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_THAT(p2.force, AnyOf(Eq(vec3(3,4,5)), Eq(-vec3(3,4,5))));
	EXPECT_EQ(p1.force, p2.force);
}

TYPED_TEST(LinkedCellsTest, TwoParticles_ConstantTypeForce_NoNeighbouringCell) {
	Environment e(forces<ConstantForce>);
	e.set_extent({2,1,0.5});
	e.set_origin({0,0,0});
	e.add_particle(make_particle(7, {0.25,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(7, {1.25,0,0}, {}, 1, ParticleState::ALIVE, 1));
	e.add_force(ConstantForce(3,4,5), to_type(7));

	auto sys = build_system(e, TypeParam::create_container(0.5), TypeParam::create_exec());
	sys.update_forces();
    auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 2u);

	auto & p1 = out[0].mass == 1 ? out[0] : out[1];
	auto & p2 = out[0].mass == 2 ? out[0] : out[1];

	// particles should not interact because they are not in neighboring cells
	EXPECT_EQ(p1.force, vec3(0,0,0));
	EXPECT_EQ(p2.force, vec3(0,0,0));
}

TYPED_TEST(LinkedCellsTest, TwoParticles_IdSpecificForce) {
    Environment e(forces<NoForce, ConstantForce>);
	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 42));
	e.add_particle(make_particle(0, {0,1,0}, {}, 1, ParticleState::ALIVE, 99));
	e.add_force(NoForce(), to_type(0));
	e.add_force(ConstantForce(-1,2,-3), between_ids(42, 99));
	e.auto_domain(2);

	auto sys = build_system(e, TypeParam::create_container(), TypeParam::create_exec());
	sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), 2u);

	EXPECT_THAT(
		out[0].force,
		AnyOf(Eq(vec3(-1,2,-3)), Eq(-vec3(-1,2,-3)))
	);
}

TYPED_TEST(LinkedCellsTest, TwoParticles_InverseSquare) {
    Environment e(forces<NoForce, Gravity>);

    e.set_extent({10,10,10});

	e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
	e.add_particle(make_particle(1, {2,0,0}, {}, 2, ParticleState::ALIVE, 1));


	e.add_force(NoForce(), to_type(0));
	e.add_force(NoForce(), to_type(1));

	e.add_force(Gravity(5.0), between_types(0, 1));

	auto sys = build_system(e, TypeParam::create_container(), TypeParam::create_exec());
	sys.update_forces();

    auto const& out = export_particles(sys);
    // find each
    const auto& pa = (out[0].mass == 1 ? out[0] : out[1]);
    const auto& pb = (out[1].mass == 2 ? out[1] : out[0]);
    // magnitude = pre * m1*m2 / r^3 = 5*1*2/(2^3)=10/8=1.25  direction from pa->pb = (2,0,0)
    // force on pa = 1.25*(2,0,0) = (2.5,0,0); on pb = (-2.5,0,0)
    EXPECT_NEAR(pa.force.x, 2.5, 1e-12);
    EXPECT_NEAR(pb.force.x, -2.5, 1e-12);
    EXPECT_EQ(pa.force.y, 0.0);
    EXPECT_EQ(pb.force.y, 0.0);
}


TYPED_TEST(LinkedCellsTest, OrbitTest) {
	constexpr double G = 1;
	constexpr double R = 1;
	constexpr double M = 1.0;
	constexpr double m = 1e-10;
	constexpr double v = G * M / R;
	constexpr double T = 2 * 3.14159265359 * v / R;

	Environment env (forces<Gravity>);

	env.add_particle(make_particle(0, {0,R,0}, {v, 0, 0}, m));
	env.add_particle(make_particle(0, {0,0,0}, {0, 0, 0}, M));

	env.add_force(Gravity(G), to_type(0));

	env.set_origin({-1.5*v,-1.5*v,0});
	env.set_extent({3*v,3*v,1});

	auto sys = build_system(env, TypeParam::create_container(v), TypeParam::create_exec());
	sys.update_forces();

	VelocityVerlet integrator(sys, monitors<OrbitMonitor>);
	integrator.add_monitor(OrbitMonitor(v, R));
	integrator.run_for_duration(0.001, T);

	auto particles = export_particles(sys);

	auto p1 =  particles[0].mass == m ? particles[0] : particles[1];
	auto p2 =  particles[0].mass == M ? particles[0] : particles[1];

	EXPECT_NEAR(p1.velocity.norm(), v, 1e-3);

	EXPECT_NEAR(p1.position.x, 0, 1e-3);
	EXPECT_NEAR(p1.position.y, R, 1e-3);
	EXPECT_EQ(p1.position.z, 0);

	EXPECT_NEAR(p1.velocity.x, v, 1e-3);
	EXPECT_NEAR(p1.velocity.y, 0, 1e-3);
	EXPECT_EQ(p1.velocity.z, 0);

	EXPECT_NEAR(p2.position.x, 0, 1e-3);
	EXPECT_NEAR(p2.position.y, 0, 1e-3);
	EXPECT_NEAR(p2.position.z, 0, 1e-3);

	EXPECT_NEAR(p2.velocity.x, 0, 1e-3);
	EXPECT_NEAR(p2.velocity.y, 0, 1e-3);
	EXPECT_NEAR(p2.velocity.z, 0, 1e-3);
}


TYPED_TEST(LinkedCellsTest, CollectIndicesInRegion) {
	// Create a simple 3x3x3 grid of particles (27 total)
	auto cuboid = ParticleCuboid{}
		.at(vec3(0.25))
		.velocity({0, 0, 0})
		.count({3, 3, 3})
		.mass(1.0)
		.spacing(1)
		.type(0);

    // Loop through different cell size hints to verify consistency
    for (double cell_size : {0.5, 1.0, 2.0, 5.0}) {

        Environment e(forces<NoForce>);
        e.set_origin({0, 0, 0});
        e.set_extent({5, 5, 5});
		e.add_particles(cuboid);
        e.add_force(NoForce(), to_type(0));

        auto sys = build_system(e, TypeParam::create_container(cell_size), TypeParam::create_exec());

        // Case 1: small inner region (should include one particle)
        {
            Domain region({0.1, 0.1, 0.1}, {0.9, 0.9, 0.9});
            auto indices = sys.query_region(region);
            ASSERT_EQ(indices.size(), 1u);
	        auto p = export_particles(sys)[indices[0]];
            EXPECT_EQ(p.position.x, 0.25);
            EXPECT_EQ(p.position.y, 0.25);
            EXPECT_EQ(p.position.z, 0.25);
        }

        // Case 2: mid region (should include all 27)
        {
            Domain region({0, 0, 0}, {5, 5, 5});
            auto indices = sys.query_region(region);
            EXPECT_EQ(indices.size(), 27u);
        }

        // Case 3: partially overlapping region
        {
            Domain region({1.5, 1.5, 1.5}, {4.5, 4.5, 4.5});
            std::vector indices = sys.query_region(region);
            EXPECT_GT(indices.size(), 0u);
            EXPECT_LT(indices.size(), 27u);

        	std::unordered_set inside (indices.begin(), indices.end());

			for (size_t id = sys.min_id(); id < sys.max_id(); id++) {
				auto p = get_particle(sys, id);
        		bool in_region = (p.position.x >= 1.5 && p.position.x <= 4.5) &&
					(p.position.y >= 1.5 && p.position.y <= 4.5) &&
					(p.position.z >= 1.5 && p.position.z <= 4.5);

        		if (inside.contains(id)) {
        			EXPECT_TRUE(in_region);
        		} else {
        			EXPECT_FALSE(in_region);
        		}
        	}
        }

        // Case 4: region completely outside
        {
            Domain region({10, 10, 10}, {12, 12, 12});
            auto indices = sys.query_region(region);
            EXPECT_TRUE(indices.empty());
        }
    }
}

// does nothing except signaling the container to be periodic
struct DummyPeriodicBoundary final : boundary::Boundary {
	static constexpr ParticleField fields = ParticleField::none;

	DummyPeriodicBoundary()
	: Boundary(0.0, false, true, false ) {}

	template<ParticleField M, particle::IsParticleAttributes U>
		void apply(auto, const core::Box &, const DomainFace) const noexcept{
	}
};

TYPED_TEST(LinkedCellsTest, PeriodicForceWrap_X) {
	// Iterate over several cell sizes (smaller, medium, larger than extent/2)
	for (double cell_size_hint : {9.9}) {
		Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
		e.set_origin({0,0,0});
		e.set_extent({10,10,10}); // domain box 10x10x10

		// Two particles, near opposite faces along x
		e.add_particle(make_particle(0, {0.5, 5, 5}, {}, 1, ParticleState::ALIVE, 0));
		e.add_particle(make_particle(0, {9.5, 5, 5}, {}, 1, ParticleState::ALIVE, 1));

		// Simple harmonic force
		e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

		// Enable periodic boundaries on both x faces
		e.set_boundaries(DummyPeriodicBoundary(), {DomainFace::XMinus, DomainFace::XPlus});

		BuildInfo mapping;
		auto sys = build_system(e, TypeParam::create_container(cell_size_hint), TypeParam::create_exec(), &mapping);
		sys.update_forces();

		auto const& out = export_particles(sys);
		ASSERT_EQ(out.size(), 2u);

		auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
		auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

		EXPECT_NEAR(p1.force.x, -1.0, 1e-12);
		EXPECT_NEAR(p2.force.x, 1.0, 1e-12);
	}
}


TYPED_TEST(LinkedCellsTest, PeriodicForceWrap_AllAxes) {
	for (double cell_size_hint : {1.0, 3.3, 9.9}) {

		Environment e(forces<Harmonic>, boundaries<DummyPeriodicBoundary>);
		e.set_origin({0,0,0});
		e.set_extent({10,10,10});

		// Particles at opposite corners
		e.add_particle(make_particle(0, {0.5, 0.5, 0.5}, {}, 1, ParticleState::ALIVE, 0));
		e.add_particle(make_particle(0, {9.5, 9.5, 9.5}, {}, 1, ParticleState::ALIVE, 1));

		e.add_force(Harmonic(1.0, 0.0, 2.0), to_type(0));

		// Enable full periodicity on all faces
		e.set_boundaries(DummyPeriodicBoundary(), {
			DomainFace::XMinus, DomainFace::XPlus,
			DomainFace::YMinus, DomainFace::YPlus,
			DomainFace::ZMinus, DomainFace::ZPlus
		});

		BuildInfo mapping;
		auto sys = build_system(e, TypeParam::create_container(cell_size_hint), TypeParam::create_exec(), &mapping);
		sys.update_forces();

		auto const& out = export_particles(sys);
		ASSERT_EQ(out.size(), 2u);

		auto p1 = get_particle_by_id(sys, mapping.id_map[0]);
		auto p2 = get_particle_by_id(sys, mapping.id_map[1]);

		// Forces must be equal and opposite
		EXPECT_EQ(p1.force, -p2.force);

		EXPECT_NEAR(p1.force.x, -1.0, 1e-12);
		EXPECT_NEAR(p1.force.y, -1.0, 1e-12);
		EXPECT_NEAR(p1.force.z, -1.0, 1e-12);

		EXPECT_NEAR(p2.force.x, 1.0, 1e-12);
		EXPECT_NEAR(p2.force.y, 1.0, 1e-12);
		EXPECT_NEAR(p2.force.z, 1.0, 1e-12);
	}
}

TYPED_TEST(LinkedCellsTest, Asymmetric_ChunkBoundaries_Counting) {
    constexpr size_t n_type0 = 20;
    constexpr size_t n_type1 = 12;

    Environment e(forces<Harmonic, NoForce>);
    e.set_extent({10, 10, 10});

    for (ParticleID i = 0; i < n_type0; ++i) {
        e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, i));
    }
    for (ParticleID i = 0; i < n_type1; ++i) {
        e.add_particle(make_particle(1, {1,0,0}, {}, 1, ParticleState::ALIVE, 100 + i));
    }

    // Spring k=1, r0=0, cutoff=1.5 (covers dist=1.0)
    e.add_force(Harmonic(1, 0, 1.5), between_types(0, 1));

    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));

    BuildInfo info;
    // Set cell size >= cutoff
    auto sys = build_system(e, TypeParam::create_container(1.5), TypeParam::create_exec(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);
    ASSERT_EQ(out.size(), n_type0 + n_type1);

    // P0 pulls P1 (+X), P1 pulls P0 (-X)
    const vec3 expected_f0 = vec3(1, 0, 0) * static_cast<double>(n_type1);
    const vec3 expected_f1 = vec3(-1, 0, 0) * static_cast<double>(n_type0);

    for (const auto& p : out) {
        if (p.type == info.type_map[0]) EXPECT_EQ(p.force, expected_f0);
        else EXPECT_EQ(p.force, expected_f1);
    }
}

TYPED_TEST(LinkedCellsTest, Asymmetric_MultiChunk_Gravity_WithCutoff) {
    constexpr size_t n_a = 10;
    constexpr size_t n_b = 10;
    constexpr double cutoff = 12.0;

    Environment e(forces<Gravity, NoForce>);
    e.set_extent({100, 100, 100});

    // Line of Type 0 at y=0
    for (ParticleID i = 0; i < n_a; ++i) {
        e.add_particle(make_particle(0, {static_cast<double>(i), 0, 0}, {}, 1.0, ParticleState::ALIVE, i));
    }
    // Line of Type 1 at y=10 (dist ~10.0)
    for (ParticleID i = 0; i < n_b; ++i) {
        e.add_particle(make_particle(1, {static_cast<double>(i), 10, 0}, {}, 1.0, ParticleState::ALIVE, 100+i));
    }

    // Gravity G=1, Cutoff=12 (covers the gap of 10)
    e.add_force(Gravity(1.0, cutoff), between_types(0, 1));
    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));

    BuildInfo info;
    auto sys = build_system(e, TypeParam::create_container(cutoff), TypeParam::create_exec(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);

    for (const auto& p : out) {
        vec3 expected_force(0,0,0);
        const int target_user_type = (p.type == info.type_map[0]) ? 1 : 0;

        for (const auto& other : out) {
            if (other.type != info.type_map[target_user_type]) continue;

            vec3 r = other.position - p.position;
            const double dist = r.norm();

            // LC only calculates if within cutoff
            if (dist > cutoff) continue;

            const double mag = (1.0 * p.mass * other.mass) / (dist * dist * dist);
            expected_force += r * mag;
        }

        EXPECT_NEAR(p.force.x, expected_force.x, 1e-10);
        EXPECT_NEAR(p.force.y, expected_force.y, 1e-10);
        EXPECT_NEAR(p.force.z, expected_force.z, 1e-10);
    }
}

TYPED_TEST(LinkedCellsTest, Asymmetric_TypeChaining) {
    Environment e(forces<Harmonic, NoForce>);
    e.set_extent({10, 10, 10});

    e.add_particle(make_particle(0, {0,0,0}, {}, 1, ParticleState::ALIVE, 0));
    e.add_particle(make_particle(1, {1,0,0}, {}, 1, ParticleState::ALIVE, 1));
    e.add_particle(make_particle(2, {1,1,0}, {}, 1, ParticleState::ALIVE, 2));

    // Cutoff 1.5 covers the distance of 1.0
    e.add_force(Harmonic(100, 0, 1.5), between_types(0, 1));
    e.add_force(Harmonic(10, 0, 1.5),  between_types(1, 2));

    e.add_force(NoForce(), to_type(0));
    e.add_force(NoForce(), to_type(1));
    e.add_force(NoForce(), to_type(2));
    e.add_force(NoForce(), between_types(0, 2));

    BuildInfo info;
    auto sys = build_system(e, TypeParam::create_container(1.5), TypeParam::create_exec(), &info);
    sys.update_forces();

    auto const& out = export_particles(sys);

    auto get_force = [&](const size_t id) {
        for (const auto& p : out) if (p.id == info.id_map[id]) return p.force;
        return vec3(0);
    };

    EXPECT_EQ(get_force(0), vec3(100, 0, 0));
    EXPECT_EQ(get_force(2), vec3(0, -10, 0));
    EXPECT_EQ(get_force(1), vec3(-100, 10, 0));
}


TYPED_TEST(LinkedCellsTest, IdBasedAccess_ReadWrite) {
	constexpr size_t N = 20;
	Environment e(forces<NoForce>);
	e.set_extent({N * 1.0, 10, 10});

	// setup: add particles
	for (ParticleID i = 0; i < N; ++i) {
		const double coord = static_cast<double>(i) + 0.5;
		e.add_particle(make_particle(0, {coord, 0.5, 0.5}, {0, 0, 0}, 1.0, ParticleState::ALIVE, i));
	}
	e.add_force(NoForce(), to_type(0));

	BuildInfo info;
	auto sys = build_system(e, TypeParam::create_container(1.0), TypeParam::create_exec(), &info);

	// test view_id (Read Only)
	for (size_t i = 0; i < N; ++i) {
		// Resolve the internal ID used by the system
		const auto sys_id = info.id_map[i];
		auto view = sys.template view_id<ParticleField::position | ParticleField::id>(sys_id);
		EXPECT_EQ(view.id, sys_id);
		EXPECT_DOUBLE_EQ(view.position.x, static_cast<double>(i) + 0.5);
	}

	// test at_id (Read/Write)
	for (size_t i = 0; i < N; ++i) {
		const auto sys_id = info.id_map[i];

		// Modify velocity using ID access
		auto ref = sys.template at_id<ParticleField::velocity>(sys_id);
		const auto val = static_cast<double>(i);
		ref.velocity = {val, val * 2, val * 3};
	}

	// test restricted_at_id (Read Verification of previous Write)
	for (size_t i = 0; i < N; ++i) {
		const auto sys_id = info.id_map[i];
		const auto val = static_cast<double>(i);

		// Verify the write persisted and can be read back via restricted interface
		auto res_view = sys.template at_id<ParticleField::velocity | ParticleField::force>(sys_id);

		EXPECT_EQ(res_view.velocity.x, val);
		EXPECT_EQ(res_view.velocity.y, val * 2);
		EXPECT_EQ(res_view.velocity.z, val * 3);
	}
}


TYPED_TEST(LinkedCellsTest, Sparse_SIMD_Mask_Check) {
	// This specifically targets the AoSoA padding logic.
	Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
	env.add_force(LennardJones(5.0, 1.0, 3.0), to_type(0));

	// Massive extent, tiny number of particles = lots of empty cells
	env.set_extent({100, 100, 100});
	env.set_origin({0, 0, 0});

	// Add just 4 particles, completely isolated from each other
	env.add_particle(make_particle(0, {10, 10, 10}, {}, 1.0, ParticleState::ALIVE, 0));
	env.add_particle(make_particle(0, {50, 50, 50}, {}, 1.0, ParticleState::ALIVE, 1));
	env.add_particle(make_particle(0, {90, 90, 90}, {}, 1.0, ParticleState::ALIVE, 2));

	// Add two particles right next to each other to trigger exactly ONE interaction
	env.add_particle(make_particle(0, {25, 25, 25}, {}, 1.0, ParticleState::ALIVE, 3));
	env.add_particle(make_particle(0, {25.5, 25, 25}, {}, 1.0, ParticleState::ALIVE, 4));

	auto sys = build_system(env, TypeParam::create_container(3.0), TypeParam::create_exec());
	sys.update_forces();

	auto const& out = export_particles(sys);
	ASSERT_EQ(out.size(), 5u);

	// Particles 0, 1, and 2 should have EXACTLY zero force.
	for (const auto& p : out) {
		if (p.id == 0 || p.id == 1 || p.id == 2) {
			EXPECT_EQ(p.force, vec3(0,0,0)) << "SIMD Padding Leak detected on isolated particle " << p.id;
		}
	}
}



TYPED_TEST(LinkedCellsTest, LinkedCells_vs_DirectSum_Parity_Open) {

	// Setup Environment with OpenBoundary
	// This ensures both LC and DS only see the "real" distances.
	Environment env(forces<LennardJones>, boundaries<OpenBoundary>);
	env.add_force(LennardJones(5.0, 1.0, 3.0), to_type(0));

	// Safe Jittered Grid Initialization
	const int nx = 10, ny = 10, nz = 10;
	const double spacing = 1.15; // Slightly larger than sigma (1.0) to stay out of the wall
	const double jitter_mag = 0.05; // Small nudge

	std::mt19937 gen(42);
	std::uniform_real_distribution<double> jitter(-jitter_mag, jitter_mag);

	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				ParticleID user_id = k * (nx * ny) + j * nx + i;

				vec3 pos = {
					0.5 + i * spacing + jitter(gen),
					0.5 + j * spacing + jitter(gen),
					0.5 + k * spacing + jitter(gen)
				};

				env.add_particle(make_particle(0, pos, {0,0,0}, 1.0, ParticleState::ALIVE, user_id));
			}
		}
	}

	CustomExecConfig<ParallelPolicy::Serial> serial_exec;

	// Build systems
	BuildInfo ds_info, lc_info;
	auto ds_system = build_system(env, DirectSum<Layout::AoS>{}, serial_exec, &ds_info);
	auto lc_system = build_system(env, TypeParam::create_container(3.0), serial_exec, &lc_info);

	ds_system.update_forces();
	lc_system.update_forces();

	// Comparison
	constexpr double epsilon = 1e-9;
	bool found_failure = false;
	for (ParticleID user_id = 0; user_id < 1000; ++user_id) {
		const auto ds_idx = ds_info.id_map[user_id];
		const auto lc_idx = lc_info.id_map[user_id];

		auto p_ds = get_particle_by_id(ds_system, ds_idx);
		auto p_lc = get_particle_by_id(lc_system, lc_idx);

		// 1. POSITION GUARD: If these aren't the same, the mapping is broken.
		double pos_diff = (p_ds.position - p_lc.position).norm();
		if (pos_diff > epsilon) {
			std::cout << "\n[  MAPPING FAIL  ] User ID: " << user_id << "\n"
					  << "  DS Internal Idx: " << ds_idx << " | LC Internal Idx: " << lc_idx << "\n"
					  << "  DS Pos: " << p_ds.position << "\n"
					  << "  LC Pos: " << p_lc.position << "\n";
			found_failure = true;
			break;
		}

		// 2. FORCE CHECK
		if (std::abs(p_ds.force.x - p_lc.force.x) > epsilon) {
			std::cout << "\n[   FORCE FAIL   ] User ID: " << user_id << "\n"
					  << "  Position: " << p_ds.position << "\n"
					  << "  DS Force: " << p_ds.force << "\n"
					  << "  LC Force: " << p_lc.force << "\n"
					  << "  Delta:    " << (p_ds.force - p_lc.force) << "\n";
			found_failure = true;
			break;
		}
	}

	ASSERT_FALSE(found_failure) << "Parity failed. See log for details.";
}






