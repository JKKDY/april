#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include "april/april.h"
#include "utils.h"

using namespace april;

// Define the "Sinks" to store test data
struct SpySinks {
    int init_call_count = 0;
    int apply_call_count = 0;
    std::vector<size_t> steps_at_apply;
    std::vector<double> times_at_apply;
};

// It holds a raw pointer to the sinks object owned by the test.
class SpyController : public controller::Controller {
public:
    SpySinks* sinks = nullptr; // Raw pointer, not an owner

    // Constructor used in the test to pass in the sinks
    explicit SpyController(shared::Trigger trig, SpySinks* sinks_ptr)
        : Controller(std::move(trig)), sinks(sinks_ptr) {}

    // Default constructor is required for the controller::controllers<...> pack
    SpyController() : Controller(Trigger::never()) {}

    // Implement init: increment counter
    template<class S>
    void init(const SystemContext<S>&) {
        if (sinks) sinks->init_call_count++;
    }

    // Implement apply: increment counter and record step/time
    template<class S>
    void apply(SystemContext<S>& ctx) {
        if (sinks) {
            sinks->apply_call_count++;
            sinks->steps_at_apply.push_back(ctx.step());
            sinks->times_at_apply.push_back(ctx.time());
        }
    }
};


// This fixture will own the SpySinks object.
class ControllerTest : public testing::Test {
protected:
    // The test fixture owns the sinks.
    SpySinks sinks;

    // Helper function to build a minimal environment
    // ready for testing.
    auto setup_environment(const Trigger& trigger) {
        // 1. Reset the sinks at the start of each test
        sinks = {};

        // 2. Set up a minimal environment
        return Environment(forces<NoForce>, controllers<SpyController>)
            .with_particle(Particle().at({}).as_type(0).with_mass(1)) // Need one particle
            .with_force(NoForce(), to_type(0))
            .with_controller(SpyController(trigger, &sinks))
            .with_extent(1,1,1);
    }
};

// Each test now builds its own system and integrator,
// ensuring the 'system' outlives the 'integrator'.

TEST_F(ControllerTest, InitIsCalledOnce) {
    const auto env = setup_environment(Trigger::never());
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(10).run();

    // Init should be called exactly once by the integrator
    EXPECT_EQ(sinks.init_call_count, 1);
}

TEST_F(ControllerTest, TriggerNever) {
    const auto env = setup_environment(Trigger::never());
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(10).run();

    EXPECT_EQ(sinks.apply_call_count, 0);
    EXPECT_TRUE(sinks.steps_at_apply.empty());
}

TEST_F(ControllerTest, TriggerAlways) {
    const auto env = setup_environment(Trigger::always());
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run(); // Runs steps 0, 1, 2, 3, 4

    EXPECT_EQ(sinks.apply_call_count, 5);
    std::vector<size_t> expected_steps = { 0, 1, 2, 3, 4 };
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerEvery3Steps) {
    const auto env = setup_environment(Trigger::every(3));
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(10).run(); // Runs 1..10. Triggers at 0, 3, 6, 9.

    EXPECT_EQ(sinks.apply_call_count, 4);
    std::vector<size_t> expected_steps = { 0, 3, 6, 9 };
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerEvery3StepsWithOffset) {
    const auto env = setup_environment(Trigger::every(3, 1)); // offset = 1
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(10).run(); // (step + 1) % 3 == 0. Triggers at 2, 5, 8.

    EXPECT_EQ(sinks.apply_call_count, 3);
    const std::vector<size_t> expected_steps = {2, 5, 8};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerAtStep5) {
    const auto env = setup_environment(Trigger::at_step(5));
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(10).run();

    EXPECT_EQ(sinks.apply_call_count, 1);
    const std::vector<size_t> expected_steps = {5};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerAfterStep4) {
    const auto env = setup_environment(Trigger::after(4));
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(7).run(); // Triggers at 4, 5, 6

    EXPECT_EQ(sinks.apply_call_count, 3);
    const std::vector<size_t> expected_steps = {4, 5, 6};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerBetweenSteps3And5) {
    const auto env = setup_environment(Trigger::between(3, 5));
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(7).run(); // Triggers at 3, 4. (end is exclusive)

    EXPECT_EQ(sinks.apply_call_count, 2);
    const std::vector<size_t> expected_steps = {3, 4};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerAfterTime) {
    const auto env = setup_environment(Trigger::after_time(0.025));
    auto system = build_system(env, DirectSum());
    // dt=0.01. Steps: 1 (t=0.01), 2 (t=0.02), 3 (t=0.03)
    // Should trigger at step 3.
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run();

    EXPECT_EQ(sinks.apply_call_count, 2);
    const std::vector<size_t> expected_steps = {3, 4};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
    EXPECT_NEAR(sinks.times_at_apply[0], 0.03, 1e-12);
}

TEST_F(ControllerTest, TriggerPeriodically) {
    const auto env = setup_environment(Trigger::periodically(0.03));
    auto system = build_system(env, DirectSum());
    // dt=0.01. Steps/Time: 1/0.01, 2/0.02, 3/0.03, 4/0.04, 5/0.05, 6/0.06, 7/0.07
    // Fires at 1 (t=0.01), 4 (t=0.04), 7 (t=0.07)
    StoermerVerlet(system).with_dt(0.01).for_steps(8).run();

    EXPECT_EQ(sinks.apply_call_count, 3);
    const std::vector<size_t> expected_steps = {0, 3, 6};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerLogicalOr) {
    const auto trigger = Trigger::at_step(2) || Trigger::at_step(4);
    const auto env = setup_environment(trigger);
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run();

    EXPECT_EQ(sinks.apply_call_count, 2);
    const std::vector<size_t> expected_steps = {2, 4};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerLogicalAnd) {
    const auto trigger = Trigger::every(2) && Trigger::after(4);
    const auto env = setup_environment(trigger);
    auto system = build_system(env, DirectSum());
    // Steps: 0, 1, 2, 3, 4, 5, 6, 7
    // every(2): 0, 2, 4, 6
    // after(4): 3, 4, 5, 6, 7
    // AND: 4, 6
    StoermerVerlet(system).with_dt(0.01).for_steps(8).run();

    EXPECT_EQ(sinks.apply_call_count, 2);
    const std::vector<size_t> expected_steps = {4, 6};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}

TEST_F(ControllerTest, TriggerLogicalNot) {
    // Trigger at every step EXCEPT step 3
    const auto trigger = Trigger::always() && !Trigger::at_step(3);
    const auto env = setup_environment(trigger);
    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run(); // 1, 2, 4, 5

    EXPECT_EQ(sinks.apply_call_count, 4);
    const std::vector<size_t> expected_steps = {0, 1, 2, 4};
    EXPECT_EQ(sinks.steps_at_apply, expected_steps);
}


class SpyController2 : public controller::Controller {
public:
    SpySinks* sinks = nullptr;
    SpyController2(shared::Trigger trig, SpySinks* sinks_ptr)
        : Controller(std::move(trig)), sinks(sinks_ptr) {}
    SpyController2() : Controller(Trigger::never()), sinks(nullptr) {}

    template<class S>
    void init(const core::SystemContext<S>&) {
        if (sinks) sinks->init_call_count++;
    }
    template<class S>
    void apply(core::SystemContext<S>& ctx) {
        if (sinks) {
            sinks->apply_call_count++;
            sinks->steps_at_apply.push_back(ctx.step());
        }
    }
};

// A controller that uses the SystemContext to find
// and modify a particle by its ID.
class ContextSpyController : public controller::Controller {
public:
    static constexpr env::FieldMask mask = env::to_field_mask(env::Field::velocity);

    ContextSpyController(shared::Trigger trig, ParticleID id)
        : Controller(std::move(trig)), target_id(id) {}

    ContextSpyController() : Controller(Trigger::never()), target_id(0) {}

    template<class S>
    void apply(core::SystemContext<S>& ctx) {
        // Use the context to get a particle by ID and modify it
        auto p = ctx.template get_particle_by_id<mask>(target_id);
        p.velocity = {100.0, 200.0, 300.0};
    }

private:
    ParticleID target_id;
};


TEST_F(ControllerTest, MultipleSameTypeControllers) {
    // Test that two controllers of the same type can be added
    // and are triggered independently.

    // Sinks for each controller
    SpySinks sinks1;
    SpySinks sinks2;

    Environment env(
        forces<NoForce>,
        boundaries<Open>,
        controllers<SpyController>, // Only one SpyController type
        fields<>
    );
    env.with_particle(Particle().at({}).as_type(0).with_mass(1))
       .with_force(NoForce(), to_type(0))
       .with_boundaries(Open(), all_faces)
       .with_extent(1,1,1);

    // Add two separate instances of SpyController
    env.with_controller(SpyController(Trigger::at_step(2), &sinks1));
    env.with_controller(SpyController(Trigger::at_step(4), &sinks2));

    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run(); // Runs 0, 1, 2, 3, 4

    // Check sinks1
    EXPECT_EQ(sinks1.apply_call_count, 1);
    EXPECT_EQ(sinks1.steps_at_apply[0], 2);

    // Check sinks2
    EXPECT_EQ(sinks2.apply_call_count, 1);
    EXPECT_EQ(sinks2.steps_at_apply[0], 4);
}


TEST_F(ControllerTest, MultipleDifferentControllers) {
    // Test that two different controller types can be added
    // and are triggered independently.

    SpySinks sinks1;
    SpySinks sinks2;

    Environment env(
        forces<NoForce>,
        boundaries<Open>,
        controllers<SpyController, SpyController2>, // Two different types
        fields<>
    );
    env.with_particle(Particle().at({}).as_type(0).with_mass(1))
       .with_force(NoForce(), to_type(0))
       .with_boundaries(Open(), all_faces)
       .with_extent(1,1,1);

    // Add one instance of each controller
    env.with_controller(SpyController(Trigger::at_step(1), &sinks1));
    env.with_controller(SpyController2(Trigger::at_step(3), &sinks2));

    auto system = build_system(env, DirectSum());
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run(); // Runs 0, 1, 2, 3, 4

    // Check sinks1
    EXPECT_EQ(sinks1.init_call_count, 1); // Init is always called
    EXPECT_EQ(sinks1.apply_call_count, 1);
    EXPECT_EQ(sinks1.steps_at_apply[0], 1);

    // Check sinks2
    EXPECT_EQ(sinks2.init_call_count, 1);
    EXPECT_EQ(sinks2.apply_call_count, 1);
    EXPECT_EQ(sinks2.steps_at_apply[0], 3);
}

TEST_F(ControllerTest, ContextAccess_ModifiesParticles) {
    // Test that a controller can use the SystemContext
    // to find and modify a particle.

    constexpr ParticleID target_id = 0;
    const vec3 target_vel = {100.0, 200.0, 300.0};

    Environment env(
        forces<NoForce>,
        boundaries<Open>,
        controllers<ContextSpyController>, // Register our new controller
        fields<>
    );
    // Add a particle with a specific ID and zero velocity
    env.with_particle(
           Particle().at({}).as_type(0).with_mass(1).with_id(target_id)
       )
       .with_force(NoForce(), to_type(0))
       .with_boundaries(Open(), all_faces)
       .with_extent(1,1,1);

    // Add the controller, set to trigger at step 2
    env.with_controller(ContextSpyController(Trigger::at_step(2), target_id));

    auto system = build_system(env, DirectSum());

    // Check initial state
    auto p1 = export_particles(system)[0];
    EXPECT_EQ(p1.velocity.x, 0.0);

    // Run the simulation
    StoermerVerlet(system).with_dt(0.01).for_steps(5).run(); // Runs 0, 1, 2, 3, 4

    // Check final state (particle was modified at step 2)
    auto p2 = export_particles(system)[0];
    EXPECT_EQ(p2.velocity.x, target_vel.x);
    EXPECT_EQ(p2.velocity.y, target_vel.y);
    EXPECT_EQ(p2.velocity.z, target_vel.z);
}