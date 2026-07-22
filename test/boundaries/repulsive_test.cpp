#include <gtest/gtest.h>


#include "../../include/april/particle/access/scalar_access.hpp"
#include "april/boundaries/boundary.hpp"
#include "april/boundaries/boundary_table.hpp"
#include "april/boundaries/repulsive.hpp"
#include "april/containers/direct_sum.hpp"
#include "../../include/april/containers/layout/layout.hpp"
#include "april/containers/linked_cells.hpp"

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

inline particle::ParticleRecord<NoParticleAttributes> make_particle(const vec3& pos) {
	particle::ParticleRecord<NoParticleAttributes> p;
	p.id = 0;
	p.position = pos;
	p.force = {0,0,0};
	p.velocity = {0,0,0};
	p.mass = 1.0;
	p.state = ParticleState::ALIVE;
	return p;
}

template<ParticleField Mask, typename RecordT>
auto make_source(RecordT& record) {
	// Determine constness based on RecordT (allows making const sources from const records)
	using UserDataT = RecordT::particle_attributes_t;

	particle::internal::ParticleSource<Mask, Mask, UserDataT> src;

	if constexpr (particle::internal::has_field_v<Mask, ParticleField::position>)     src.position     = &record.position;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::velocity>)     src.velocity     = &record.velocity;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::force>)        src.force        = &record.force;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::old_position>) src.old_position = &record.old_position;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::mass>)         src.mass         = &record.mass;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::state>)        src.state        = &record.state;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::type>)         src.type         = &record.type;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::id>)           src.id           = &record.id;
	if constexpr (particle::internal::has_field_v<Mask, ParticleField::attributes>)    src.attributes    = &record.attributes;

	return src;
}

// Direct application test should add constant force
TEST(RepulsiveBoundaryTest, Apply_AddsInwardForce) {
	ConstantForce f{5.0, 10.0};
	const RepulsiveBoundary rep(f);
	constexpr ParticleField Mask = RepulsiveBoundary<ConstantForce>::fields;

	auto p = make_particle({9.5,5,5});
	auto src = make_source<Mask>(p);
	particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);

	const core::Box box({0,0,0}, {10,10,10});

	rep.apply(ref, box, DomainFace::XPlus);

	EXPECT_NEAR(p.force.x, -5.0, 1e-12);
	EXPECT_NEAR(p.force.y,  0.0, 1e-12);
	EXPECT_NEAR(p.force.z,  0.0, 1e-12);

	// X- face should push opposite direction
	p.force = {0,0,0};
	rep.apply(ref, box, DomainFace::XMinus);
	EXPECT_NEAR(p.force.x, +5.0, 1e-12);
}

// Topology sanity: inside region, not coupled, no force wrap, np position change
TEST(RepulsiveBoundaryTest, Topology_IsInsideAndNoChangesPosition) {
	ConstantForce f{1.0, 3.0};
	const RepulsiveBoundary rep(f);

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
	std::variant<RepulsiveBoundary<ConstantForce>> variant = RepulsiveBoundary(f);
	Domain domain({0,0,0}, {10,10,10});

	constexpr ParticleField Mask = RepulsiveBoundary<ConstantForce>::fields;

	auto compiled = boundary::internal::compile_boundary(variant, core::Box::from_domain(domain), DomainFace::YMinus);

	auto p = make_particle({5,0.3,5});
	auto src = make_source<Mask>(p);
	particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);

	const core::Box box({0,0,0}, {10,10,10});

	compiled.dispatch([&](auto && bc) {
		bc.apply(ref, box, DomainFace::YMinus);
	});

	EXPECT_NEAR(p.force.y, +2.0, 1e-12)
		<< "Force on Y- face should push inward (+Y direction).";
	EXPECT_NEAR(p.force.x, 0.0, 1e-12);
	EXPECT_NEAR(p.force.z, 0.0, 1e-12);
}


// System-level pipeline test
template <class ContainerT>
class RepulsiveBoundarySystemTestT : public testing::Test {};

using ContainerTypes = testing::Types<
    DirectSum<Layout::AoS>, DirectSum<Layout::SoA>, DirectSum<Layout::AoSoA<>>,
    LinkedCells<Layout::SoA>, LinkedCells<Layout::SoA>, LinkedCells<Layout::AoSoA<>>
>;
TYPED_TEST_SUITE(RepulsiveBoundarySystemTestT, ContainerTypes);


TYPED_TEST(RepulsiveBoundarySystemTestT, EachFace_AppliesInwardForce) {
	ConstantForce f{3.0, 5.0};
	Environment env(forces<NoForce>, boundaries<RepulsiveBoundary<ConstantForce>>);
	env.set_origin({0,0,0});
	env.set_extent({10,10,10});
	env.add_interaction(NoForce{}, to_type(0));

	// Particles near each face
	env.add_particle(make_particle(0, {0.5,5,5}, {}, 1, ParticleState::ALIVE, 0)); // X-
	env.add_particle(make_particle(0, {9.5,5,5}, {}, 1, ParticleState::ALIVE, 1)); // X+
	env.add_particle(make_particle(0, {5,0.5,5}, {}, 1, ParticleState::ALIVE, 2)); // Y-
	env.add_particle(make_particle(0, {5,9.5,5}, {}, 1, ParticleState::ALIVE, 3)); // Y+
	env.add_particle(make_particle(0, {5,5,0.5}, {}, 1, ParticleState::ALIVE, 4)); // Z-
	env.add_particle(make_particle(0, {5,5,9.5}, {}, 1, ParticleState::ALIVE, 5)); // Z+

	env.set_boundaries(std::array{
		RepulsiveBoundary(f), RepulsiveBoundary(f),
		RepulsiveBoundary(f), RepulsiveBoundary(f),
		RepulsiveBoundary(f), RepulsiveBoundary(f)
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
    WallForce::ExponentialForce exp_force{10.0, 2.0, 10.0};
    const RepulsiveBoundary rep(exp_force);
    constexpr ParticleField Mask = RepulsiveBoundary<WallForce::ExponentialForce>::fields;

    // Particle 1.0 unit away from 0.0 (XMinus wall)
    auto p = make_particle({1.0, 5, 5});
    auto src = make_source<Mask>(p);
    particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);
    const core::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, DomainFace::XMinus);

    // Expected: 10 * exp(-0.5) = 10 * 0.60653... = 6.0653...
    double expected = 10.0 * std::exp(-0.5);

    // Direction: XMinus wall pushes INWARD (Positive X)
    EXPECT_NEAR(p.force.x, expected, 1e-6);
}

TEST(RepulsiveBoundaryTest, PowerLawForce_CalculatesCorrectly) {
    // A=2.0, n=2.0, rc=10
    // Formula: 2.0 / d^2
    WallForce::PowerLawForce pow_force{2.0, 2.0, 10.0};
    const RepulsiveBoundary rep(pow_force);
    constexpr ParticleField Mask = RepulsiveBoundary<WallForce::PowerLawForce>::fields;

    // Particle 2.0 units away from 0.0
    auto p = make_particle({2.0, 5, 5});
    const auto src = make_source<Mask>(p);
    particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);
    const core::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, DomainFace::XMinus);

    // Expected: 2.0 / (2.0^2) = 0.5
    EXPECT_NEAR(p.force.x, 0.5, 1e-6);
}

TEST(RepulsiveBoundaryTest, LennardJones93Force_CalculatesCorrectly) {
    // eps=1, sigma=1, rc=5
    WallForce::LennardJones93Force lj93{1.0, 1.0, 5.0};
    const RepulsiveBoundary rep(lj93);
    constexpr ParticleField Mask = decltype(rep)::fields;

    // Distance = 1.0 (sigma)
    // Formula: 4*eps * (3*(s/r)^3 - 9*(s/r)^9)
    // At r=s: 4 * (3 - 9) = -24
    auto p = make_particle({1.0, 5, 5});
    const auto src = make_source<Mask>(p);
    particle::internal::ScalarParticleRef<Mask,Mask, NoParticleAttributes> ref(src);
    const core::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, DomainFace::XMinus);

    // Note: The formula provided in the snippet returns a negative value (-24) at sigma.
    // Repulsive::apply does: force += direction * magnitude.
    // XMinus direction is +1. Magnitude is -24. Result force X should be -24.
    // (This implies the 9-3 potential is attractive at this distance).
    EXPECT_NEAR(p.force.x, -24.0, 1e-6);
}

TEST(RepulsiveBoundaryTest, AdhesiveLJForce_IsAlwaysRepulsive) {
    // This force uses std::abs, so it should always push away from the wall
    WallForce::AdhesiveLJForce adj_lj{1.0, 1.0, 5.0};
    const RepulsiveBoundary rep(adj_lj);
    constexpr ParticleField Mask = decltype(rep)::fields;

    // At sigma (1.0), standard LJ Force is 24 * eps * (2 - 1) = 24.
    auto p = make_particle({1.0, 5, 5});
    auto src = make_source<Mask>(p);
    particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);
    const core::Box box({0,0,0}, {10,10,10});

    rep.apply(ref, box, DomainFace::XMinus);

    EXPECT_GT(p.force.x, 0.0) << "AdhesiveLJ should return positive magnitude, pushing X+";
    EXPECT_NEAR(p.force.x, 24.0, 1e-6);
}

// -----------------------------------------------------------------------------
// Halo Functionality Test
// -----------------------------------------------------------------------------

TEST(RepulsiveBoundaryTest, Halo_DoublesTheDistance) {
    LinearIdentityForce lin_force; // Returns distance as magnitude
    constexpr ParticleField Mask = RepulsiveBoundary<LinearIdentityForce>::fields;
    const core::Box box({0,0,0}, {10,10,10});

    // Case 1: Halo OFF
    {
        // simulate_halo = false
        const RepulsiveBoundary rep_no_halo(lin_force, false);

        // Particle at distance 2.0 from X- wall (pos=2.0)
        auto p = make_particle({2.0, 5, 5});
        auto src = make_source<Mask>(p);
        particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);

        rep_no_halo.apply(ref, box, DomainFace::XMinus);

        // Force magnitude = distance = 2.0. Direction XMinus is +1.
        EXPECT_NEAR(p.force.x, 2.0, 1e-12);
    }

    // Case 2: Halo ON
    {
        // simulate_halo = true
        const RepulsiveBoundary rep_halo(lin_force, true);

        // Particle at distance 2.0 from X- wall (pos=2.0)
        auto p = make_particle({2.0, 5, 5});
        auto src = make_source<Mask>(p);
        particle::internal::ScalarParticleRef<Mask, Mask, NoParticleAttributes> ref(src);

        rep_halo.apply(ref, box, DomainFace::XMinus);

        // Internal distance becomes distance * 2 = 4.0.
        // Force magnitude = 4.0. Direction XMinus is +1.
        EXPECT_NEAR(p.force.x, 4.0, 1e-12) << "Halo should double the effective distance passed to the force";
    }
}

