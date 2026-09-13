#include <april/april.hpp>

using namespace april;

int main() {
    // Particle types are arbitrary integer labels used to select
    // which interactions apply.
    constexpr int DEFAULT = 0;

    auto sun = Particle()
        .at(0, 0, 0)
        .with_mass(1.0)
        .as_type(DEFAULT);

    auto planet = Particle()
        .at(1, 0, 0)
        .with_velocity(0, 1, 0)
        .with_mass(1e-3)
        .as_type(DEFAULT);

    auto moon = Particle()
        .at(1.05, 0, 0)
        .with_velocity(0, 1.2, 0)
        .with_mass(1e-6)
        .as_type(DEFAULT);

    auto env = Environment(
            interactions<Gravity>,
            boundaries<OpenBoundary>
        )
        .with_particles({sun, planet, moon})
        .with_interaction(Gravity(), to_type(DEFAULT))
        .with_boundaries(OpenBoundary(), all_faces);

    auto container = DirectSum();
    auto system = build_system(env, container);

    VelocityVerlet(system, monitors<BinaryOutput>)
        .with_monitor(BinaryOutput(Trigger::every(40), "output/"))
        .with_dt(0.005)
        .for_duration(200)
        .run();
}