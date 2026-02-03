#include <gtest/gtest.h>


#include "april/particle/access.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/boundaries/boundary_table.hpp"
#include "april/boundaries/repulsive.hpp"
#include "utils.h"

using namespace april;

struct ConstantForce {
	double value;
	double rc;
	explicit ConstantForce(const double value = 1.0, const double rc = 10.0)
		: value(value), rc(rc) {}

	[[nodiscard]] double cutoff() const noexcept { return rc; }

	[[nodiscard]] double apply(const double dist) const noexcept {
		return (dist <= rc) ? value : 0.0;
	}
};

// A "Spy" force that simply returns the distance passed to it.
// Used to verify if simulate_halo doubles the distance correctly.
struct LinearIdentityForce {
	double rc;
	explicit LinearIdentityForce(double rc = 100.0) : rc(rc) {}

	[[nodiscard]] double cutoff() const noexcept { return rc; }

	[[nodiscard]] double apply(const double dist) const noexcept {
		return dist;
	}
};

inline env::internal::ParticleRecord<env::NoUserData> make_particle(const vec3& pos) {
	env::internal::ParticleRecord<env::NoUserData> p;
	p.id = 0;
	p.position = pos;
	p.force = {0,0,0};
	p.velocity = {0,0,0};
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

template<env::FieldMask Mask, typename RecordT>
auto make_source(RecordT& record) {
	// Determine constness based on RecordT (allows making const sources from const records)
	constexpr bool IsConst = std::is_const_v<RecordT>;
	using UserDataT = typename RecordT::user_data_t;

	env::ParticleSource<Mask, UserDataT, IsConst> src;

	if constexpr (env::has_field_v<Mask, env::Field::position>)     src.position     = &record.position;
	if constexpr (env::has_field_v<Mask, env::Field::velocity>)     src.velocity     = &record.velocity;
	if constexpr (env::has_field_v<Mask, env::Field::force>)        src.force        = &record.force;
	if constexpr (env::has_field_v<Mask, env::Field::old_position>) src.old_position = &record.old_position;
	if constexpr (env::has_field_v<Mask, env::Field::mass>)         src.mass         = &record.mass;
	if constexpr (env::has_field_v<Mask, env::Field::state>)        src.state        = &record.state;
	if constexpr (env::has_field_v<Mask, env::Field::type>)         src.type         = &record.type;
	if constexpr (env::has_field_v<Mask, env::Field::id>)           src.id           = &record.id;
	if constexpr (env::has_field_v<Mask, env::Field::user_data>)    src.user_data    = &record.user_data;

	return src;
}

// Direct application test should add constant force
TEST(RepulsiveBoundaryTest, Apply_AddsInwardForce) {
	ConstantForce f{5.0, 10.0};
	const Repulsive rep(f);
	constexpr env::FieldMask Mask = Repulsive<ConstantForce>::fields;

	auto p = make_particle({9.5,5,5});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	const env::Box box({0,0,0}, {10,10,10});

	rep.apply(ref, box, Face::XPlus);

	EXPECT_NEAR(p.force.x, -5.0, 1e-12);
	EXPECT_NEAR(p.force.y,  0.0, 1e-12);
	EXPECT_NEAR(p.force.z,  0.0, 1e-12);

	// X- face should push opposite direction
	p.force = {0,0,0};
	rep.apply(ref, box, Face::XMinus);
	EXPECT_NEAR(p.force.x, +5.0, 1e-12);
}

// Topology sanity: inside region, not coupled, no force wrap, np position change
TEST(RepulsiveBoundaryTest, Topology_IsInsideAndNoChangesPosition) {
	ConstantForce f{1.0, 3.0};
	const Repulsive rep(f);

	const auto& topology = rep.topology;
	EXPECT_GT(topology.boundary_thickness, 0.0)
		<< "Repulsive boundaries operate inside the domain (positive thickness).";
	EXPECT_FALSE(topology.couples_axis);
	EXPECT_FALSE(topology.force_wrap);
	EXPECT_FALSE(topology.may_change_particle_position);
}

// CompiledBoundary apply test
TEST(RepulsiveBoundaryTest, CompiledBoundary_Apply_AddsInwardForce) {
	ConstantForce f{2.0, 5.0};
	std::variant<Repulsive<ConstantForce>> variant = Repulsive(f);
	env::Domain domain({0,0,0}, {10,10,10});

	constexpr env::FieldMask Mask = Repulsive<ConstantForce>::fields;

	auto compiled = boundary::internal::compile_boundary(variant, env::Box::from_domain(domain), Face::YMinus);

	auto p = make_particle({5,0.3,5});
	auto src = make_source<Mask>(p);
	env::ParticleRef<Mask, env::NoUserData> ref(src);

	const env::Box box({0,0,0}, {10,10,10});

	compiled.dispatch([&](auto && bc) {
		bc.apply(ref, box, Face::YMinus);
	});

	EXPECT_NEAR(p.force.y, +2.0, 1e-12)
		<< "Force on Y- face should push inward (+Y direction).";
	EXPECT_NEAR(p.force.x, 0.0, 1e-12);
	EXPECT_NEAR(p.force.z, 0.0, 1e-12);
}


// System-level pipeline test
template <class ContainerT>
class RepulsiveBoundarySystemTestT : public testing::Test {};
using ContainerTypes = testing::Types<DirectSumAoS, DirectSumSoA, LinkedCellsAoS, LinkedCellsSoA>;
TYPED_TEST_SUITE(RepulsiveBoundarySystemTestT, ContainerTypes);


TYPED_TEST(RepulsiveBoundarySystemTestT, EachFace_AppliesInwardForce) {
	ConstantForce f{3.0, 5.0};
	Environment env(forces<NoForce>, boundary::boundaries<Repulsive<ConstantForce>>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_force(NoForce{}, to_type(0));

	// Particles near each face
	env.add_particle(make_particle(0, {0.5,5,5}, {}, 1, ParticleState::ALIVE, 0)); // X-
	env.add_particle(make_particle(0, {9.5,5,5}, {}, 1, ParticleState::ALIVE, 1)); // X+
	env.add_particle(make_particle(0, {5,0.5,5}, {}, 1, ParticleState::ALIVE, 2)); // Y-
	env.add_particle(make_particle(0, {5,9.5,5}, {}, 1, ParticleState::ALIVE, 3)); // Y+
	env.add_particle(make_particle(0, {5,5,0.5}, {}, 1, ParticleState::ALIVE, 4)); // Z-
	env.add_particle(make_particle(0, {5,5,9.5}, {}, 1, ParticleState::ALIVE, 5)); // Z+

	env.set_boundaries(std::array{
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f),
		Repulsive(f), Repulsive(f)
	});

	BuildInfo mappings;
	auto sys = build_system(env, TypeParam(), &mappings);

	// Apply boundaries
	sys.rebuild_structure();
	sys.apply_boundary_conditions();

	// Expected inward forces per face
	const std::array expected = {
		vec3{+3,0,0}, vec3{-3,0,0},
		vec3{0,+3,0}, vec3{0,-3,0},
		vec3{0,0,+3}, vec3{0,0,-3}
	};

	for (int uid = 0; uid < 6; ++uid) {
		auto iid = mappings.id_map.at(uid);
		const auto p = get_particle_by_id(sys, iid);
		EXPECT_NEAR(p.force.x, expected[uid].x, 1e-12);
		EXPECT_NEAR(p.force.y, expected[uid].y, 1e-12);
		EXPECT_NEAR(p.force.z, expected[uid].z, 1e-12);
	}
}


TEST(RepulsiveBoundaryTest, ExponentialForce_CalculatesCorrectly) {
    // A=10, lambda=2.0, rc=10
    // Formula: 10 * exp(-d / 2.0)
    boundary::ExponentialForce exp_force{10.0, 2.0, 10.0};
    const boundary::Repulsive rep(exp_force);
    constexpr env::FieldMask Mask = boundary::Repulsive<boundary::ExponentialForce>::fields;

    // Particle 1.0 unit away from 0.0 (XMinus wall)
    auto p = make_particle({1.0, 5, 5});
    auto src = make_source<Mask>(p);
    env::ParticleRef<Mask, env::NoUserData> ref(src);
    const env::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, Face::XMinus);

    // Expected: 10 * exp(-0.5) = 10 * 0.60653... = 6.0653...
    double expected = 10.0 * std::exp(-0.5);

    // Direction: XMinus wall pushes INWARD (Positive X)
    EXPECT_NEAR(p.force.x, expected, 1e-6);
}

TEST(RepulsiveBoundaryTest, PowerLawForce_CalculatesCorrectly) {
    // A=2.0, n=2.0, rc=10
    // Formula: 2.0 / d^2
    boundary::PowerLawForce pow_force{2.0, 2.0, 10.0};
    const boundary::Repulsive rep(pow_force);
    constexpr env::FieldMask Mask = boundary::Repulsive<boundary::PowerLawForce>::fields;

    // Particle 2.0 units away from 0.0
    auto p = make_particle({2.0, 5, 5});
    auto src = make_source<Mask>(p);
    env::ParticleRef<Mask, env::NoUserData> ref(src);
    const env::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, Face::XMinus);

    // Expected: 2.0 / (2.0^2) = 0.5
    EXPECT_NEAR(p.force.x, 0.5, 1e-6);
}

TEST(RepulsiveBoundaryTest, LennardJones93Force_CalculatesCorrectly) {
    // eps=1, sigma=1, rc=5
    boundary::LennardJones93Force lj93{1.0, 1.0, 5.0};
    const boundary::Repulsive rep(lj93);
    constexpr env::FieldMask Mask = decltype(rep)::fields;

    // Distance = 1.0 (sigma)
    // Formula: 4*eps * (3*(s/r)^3 - 9*(s/r)^9)
    // At r=s: 4 * (3 - 9) = -24
    auto p = make_particle({1.0, 5, 5});
    auto src = make_source<Mask>(p);
    env::ParticleRef<Mask, env::NoUserData> ref(src);
    const env::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, Face::XMinus);

    // Note: The formula provided in the snippet returns a negative value (-24) at sigma.
    // Repulsive::apply does: force += direction * magnitude.
    // XMinus direction is +1. Magnitude is -24. Result force X should be -24.
    // (This implies the 9-3 potential is attractive at this distance).
    EXPECT_NEAR(p.force.x, -24.0, 1e-6);
}

TEST(RepulsiveBoundaryTest, AdhesiveLJForce_IsAlwaysRepulsive) {
    // This force uses std::abs, so it should always push away from the wall
    boundary::AdhesiveLJForce adj_lj{1.0, 1.0, 5.0};
    const boundary::Repulsive rep(adj_lj);
    constexpr env::FieldMask Mask = decltype(rep)::fields;

    // At sigma (1.0), standard LJ Force is 24 * eps * (2 - 1) = 24.
    auto p = make_particle({1.0, 5, 5});
    auto src = make_source<Mask>(p);
    env::ParticleRef<Mask, env::NoUserData> ref(src);
    const env::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, Face::XMinus);

    EXPECT_GT(p.force.x, 0.0) << "AdhesiveLJ should return positive magnitude, pushing X+";
    EXPECT_NEAR(p.force.x, 24.0, 1e-6);
}

// -----------------------------------------------------------------------------
// Halo Functionality Test
// -----------------------------------------------------------------------------

TEST(RepulsiveBoundaryTest, Halo_DoublesTheDistance) {
    LinearIdentityForce lin_force; // Returns distance as magnitude
    constexpr env::FieldMask Mask = boundary::Repulsive<LinearIdentityForce>::fields;
    const env::Box box({0,0,0}, {10,10,10});

    // Case 1: Halo OFF
    {
        // simulate_halo = false
        const boundary::Repulsive rep_no_halo(lin_force, false);

        // Particle at distance 2.0 from X- wall (pos=2.0)
        auto p = make_particle({2.0, 5, 5});
        auto src = make_source<Mask>(p);
        env::ParticleRef<Mask, env::NoUserData> ref(src);

        rep_no_halo.apply(ref, box, Face::XMinus);

        // Force magnitude = distance = 2.0. Direction XMinus is +1.
        EXPECT_NEAR(p.force.x, 2.0, 1e-12);
    }

    // Case 2: Halo ON
    {
        // simulate_halo = true
        const boundary::Repulsive rep_halo(lin_force, true);

        // Particle at distance 2.0 from X- wall (pos=2.0)
        auto p = make_particle({2.0, 5, 5});
        auto src = make_source<Mask>(p);
        env::ParticleRef<Mask, env::NoUserData> ref(src);

        rep_halo.apply(ref, box, Face::XMinus);

        // Internal distance becomes distance * 2 = 4.0.
        // Force magnitude = 4.0. Direction XMinus is +1.
        EXPECT_NEAR(p.force.x, 4.0, 1e-12) << "Halo should double the effective distance passed to the force";
    }
}