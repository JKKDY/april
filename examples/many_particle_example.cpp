#include <april/april.hpp>

using namespace april;

int main() {
    constexpr int DEFAULT = 0;

    auto blob = ParticleCuboid()
        .at(0, 0, 10)
        .count(10, 10, 10)
        .spacing(1.2)
        .mass(1.0)
        .type(DEFAULT)
        .thermal([](vec3 /*position*/) {
            constexpr double avg_vel = 1.0;
            return math::maxwell_boltzmann_velocity(avg_vel);
        });

    auto env = Environment(
            interactions<LennardJones>,
            boundaries<ReflectiveBoundary>,
            fields<UniformField>
        )
        .with_particles(blob)
        .with_extent(30, 30, 50)
        .with_interaction(
            LennardJones(3, 1),
            to_type(DEFAULT)
        )
        .with_field(
            UniformField({0.0, 0.0, -5.0})
        )
        .with_boundaries(
            ReflectiveBoundary(),
            all_faces
        );

    // Use linked cells for interaction traversal with SoA storage.
    auto container = LinkedCells<Layout::SoA>();
    auto system = build_system(env, container);

    VelocityVerlet(
            system,
            monitors<ProgressBar, BinaryOutput>
        )
        .with_monitor(
            ProgressBar(Trigger::every(50))
        )
        .with_monitor(
            BinaryOutput(Trigger::every(50), "output/")
        )
        .with_dt(0.001)
        .for_duration(10)
        .run();
}
